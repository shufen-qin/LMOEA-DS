function [Pop,FE] = ComplementaryEnvironmentalSelection(Pop,Arc,W,FE,problem,Boundary)

[N,m] = size(W);
D = size(Pop,2)-m;

%% The First Phase
OffPopX = GAonce(Pop(:,1:D),Arc(:,1:D),Boundary);
OffPopX = unique(OffPopX,'rows');

% Evaluation
obj = LSMOP('value',problem,m,OffPopX,D);
FE = FE + size(obj,1);
OffPop = [OffPopX,obj];

% Environmental Selection
PopCom = [Pop;Arc;OffPop];

[associate,Cosinemax] = Assign(PopCom(:,D+1:end),W);
Ns = length(unique(associate)');

if Ns < (2/3)*N
    Pop = EnvironmentalSelectionNSGAII(PopCom,N,m);
else
    Pop = EnvironmentSelectionDecomposition(PopCom,associate,Cosinemax,m);
end

%% The Second Phase
OffPopX = GA(Pop(:,1:D),Boundary,N);
obj = LSMOP('value',problem,m,OffPopX,D);
FE = FE + size(obj,1);
OffPop = [OffPopX,obj];
PopCom = [Pop;OffPop];

[associate2,Cosinemax2] = Assign(PopCom(:,D+1:end),W);
Ns = length(unique(associate2)');

if Ns < (2/3)*N
    Pop = EnvironmentalSelectionNSGAII(PopCom,N,m);
else
    Pop = EnvironmentSelectionDecomposition(PopCom,associate2,Cosinemax2,m);
end

end

function [associate,Cosinemax] = Assign(Obj,RefV)
    % Assign individuals to the reference vectors
    Obj = (Obj-repmat(min(Obj),size(Obj,1),1))./(repmat(max(Obj),size(Obj,1),1)-repmat(min(Obj),size(Obj,1),1));
    Cosine = 1-pdist2(Obj,RefV,'cosine');
    [Cosinemax,associate] = max(Cosine,[],2);
end
