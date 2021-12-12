function Offspring = GAonce(MatingPool,PopS,Boundary)
% This function is mainly used in the first reproduction.

[N,D] = size(MatingPool);
RandList = randperm(N);
MatingPool = MatingPool(RandList, :);

DisC = 20;
DisM = 20;
Ns = size(PopS,1);
Offspring=[];

for i = 1: size(MatingPool,1)
    
    k = randi([1,Ns],1);
    Pop = PopS(k,:);
    Parent = [MatingPool(i,:);Pop];
    Offspringtempt = GAhalf(Parent,Boundary,DisC,DisM,'GA');
    Offspring = [Offspring;Offspringtempt];
end
end

function Offspring = GAhalf(Parent,Boundary,DisC,DisM,type)
% This function includes the SBX crossover operator and the polynomial
% mutatoion operator.
N = size(Parent,1);
D = size(Boundary,2);
ProC = 0.9;
ProM = 1/D;

Offspring = zeros(N-1,D);

MaxParent = size(Parent,1);

MaxValue = Boundary(1,:);
MinValue = Boundary(2,:);

switch type
    case 'GA'
        % simularity binary crossover
        beta = zeros(1,D);
        miu  = rand(1,D);
        beta(miu<=0.5) = (2*miu(miu<=0.5)).^(1/(DisC+1));
        beta(miu>0.5)  = (2-2*miu(miu>0.5)).^(-1/(DisC+1));
        beta = beta.*(-1).^randi([0,1],1,D);
        beta(rand(1,D)>ProC) = 1;
        if rand <0.5
            Offspring(1,:)   = (Parent(1,:)+Parent(2,:))/2 ...
                +beta.*(Parent(1,:)-Parent(2,:))/2;
        else
            Offspring(1,:) = (Parent(1,:)+Parent(2,:))/2 ...
                -beta.*(Parent(1,:)-Parent(2,:))/2;
        end
        % polynomial mutation
        k    = rand(MaxParent-1,D);
        miu  = rand(MaxParent-1,D);
        Temp = k<=ProM & miu<0.5;
        Offspring(Temp) = Offspring(Temp)+(MaxValue(Temp)-MinValue(Temp)).* ...
            ((2.*miu(Temp)+(1-2.*miu(Temp)).*(1-(Offspring(Temp)-MinValue(Temp))./ ...
            (MaxValue(Temp)-MinValue(Temp))).^(DisM+1)).^(1/(DisM+1))-1);
        Temp = k<=ProM & miu>=0.5;
        Offspring(Temp) = Offspring(Temp)+(MaxValue(Temp)-MinValue(Temp)).* ...
            (1-(2.*(1-miu(Temp))+2.*(miu(Temp)-0.5).*(1-(MaxValue(Temp)-Offspring(Temp))./ ...
            (MaxValue(Temp)-MinValue(Temp))).^(DisM+1)).^(1/(DisM+1)));
    case 'SBX'
        % simularity binary crossover
        beta = zeros(1,D);
        miu  = rand(1,D);
        beta(miu<=0.5) = (2*miu(miu<=0.5)).^(1/(DisC+1));
        beta(miu>0.5)  = (2-2*miu(miu>0.5)).^(-1/(DisC+1));
        beta = beta.*(-1).^randi([0,1],1,D);
        beta(rand(1,D)>ProC) = 1;
        Offspring(1,:)   = (Parent(1,:)+Parent(2,:))/2+beta.*(Parent(1,:)-Parent(2,:))/2;
        Offspring(2,:) = (Parent(1,:)+Parent(2,:))/2-beta.*(Parent(1,:)-Parent(2,:))/2;  
    case 'PM'
        % polynomial mutation 
        k    = rand(MaxParent,D);
        miu  = rand(MaxParent,D);
        Temp = k<=ProM & miu<0.5;
        Parent(Temp) = Parent(Temp)+(MaxValue(Temp)-MinValue(Temp)).* ...
            ((2.*miu(Temp)+(1-2.*miu(Temp)).*(1-(Parent(Temp)-MinValue(Temp))./ ...
            (MaxValue(Temp)-MinValue(Temp))).^(DisM+1)).^(1/(DisM+1))-1);
        Temp = k<=ProM & miu>=0.5;
        Parent(Temp) = Parent(Temp)+(MaxValue(Temp)-MinValue(Temp)).* ...
            (1-(2.*(1-miu(Temp))+2.*(miu(Temp)-0.5).*(1-(MaxValue(Temp)-Parent(Temp))./ ...
            (MaxValue(Temp)-MinValue(Temp))).^(DisM+1)).^(1/(DisM+1)));
        Offspring = Parent;
end
% ±ß½çÔ¼Êø
Offspring(Offspring>MaxValue) = MaxValue(Offspring>MaxValue);
Offspring(Offspring<MinValue) = MinValue(Offspring<MinValue);
end