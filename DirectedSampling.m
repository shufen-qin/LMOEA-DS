function [Arc,FE] = DirectedSampling(Pop,Ns,W,FE,Boundary,problem)
% This a function of directed sampling based method to obtain some guiding
% solutions.

%--------------------------------------------------------------------------
[~,m] = size(W);
D = size(Pop,2)-m;  
Obj = Pop(:,D+1:D+m);
PopX = Pop(:,1:D);
Z = zeros(1,m); 

%% Define the search directions

% Classtering the reference vectors
 BW = eye(m,m);
 BW(BW==0) = 10e-7;
 
[~,Center,~,~] = kmeans(W,10);
RefW = [BW;Center];
VD = size(RefW,1);

[Best] = SelectingDirectedSolution(Obj,RefW);
BestX = PopX(Best,:); 
upper = Boundary(1,:);lower = Boundary(2,:);

Directnorm = [sqrt(sum((BestX - repmat(lower,VD,1)).^2,2));sqrt(sum((repmat(upper,VD,1) - BestX).^2,2))];
Direct = [BestX - repmat(lower,VD,1);repmat(upper,VD,1) - BestX]./repmat(Directnorm,1,D); 

%% Generate sampling solutions
vmax = sqrt(sum((Boundary(1,:)-Boundary(2,:)).^2,2));
vmin = 0;
VD = 2*VD;
v0 = vmin + rand(Ns,VD)*(vmax-vmin);
[PopNew,FE] = GeneratingSampledSolution(v0,Direct,Z,FE,Boundary,problem);
Arc = PopNew((NDSort(PopNew(:,D+1:D+m),1)==1),:); 

end

function [PopNew,FE] = GeneratingSampledSolution(v0,Direct,Z,FE,Boundary,problem)
%Generate some sampling solutions along with the guided directions
[Ns,VD] = size(v0);VD = VD/2;
D = size(Boundary,2);m = size(Z,2);
upper = Boundary(1,:);lower = Boundary(2,:);
PopNew = [];
for i = 1:Ns
    PopX = [repmat(lower,VD,1) + repmat(v0(i,1:VD)',1,D).* Direct(1:VD,:);...
        repmat(upper,VD,1) - repmat(v0(i,VD+1:end)',1,D).* Direct(VD+1:end,:)];
    PopX = max(min(repmat(upper,size(PopX,1),1),PopX),repmat(lower,size(PopX,1),1));
    
    obj = LSMOP('value',problem,m,PopX,D);
    FE = FE + size(obj,1);
    Pop = [PopX obj];
    PopNew = [PopNew; Pop];
end
end

function [Best] = SelectingDirectedSolution(Obj,W)
%Select some promising solutions to generate the guided direction.
Nr = size(W,1);
np = size(Obj,1);
Obj = (Obj-repmat(min(Obj),np,1))./(repmat(max(Obj),np,1)-repmat(min(Obj),np,1));
Best = zeros(Nr,1);

%Assign individuals 
 Cosinetemp = pdist2(Obj,W,'cosine');
 Cosine = 1-Cosinetemp;
 [~,associate] = max(Cosine,[],2);

Indflag = zeros(size(Obj,1),1);
 current = cell(Nr,1);
for i = 1:Nr
    current{i,1} = find(associate == i);
end

for i= 1:Nr
    if length(current{i,1})>1
        normf = sqrt(sum(Obj(current{i,1},:).^2,2));
        normW = sqrt(sum(repmat(W(i,:),length(current{i,1}),1).^2,2));
        Cosinefw = sum(Obj(current{i,1},:).*repmat(W(i,:),length(current{i,1}),1),2)./normW./normf;
        d1 = normf .* Cosinefw;
        [~,ind] = sort(d1,'ascend');
        Best(i,1) = current{i,1}(ind(1));
        Indflag(current{i,1}(ind(1)),1) = 1;
    elseif length(current{i,1})==1
        Best(i,1) = current{i,1}(1);
        Indflag(current{i,1}(1),1) = 1;
    end
end
for i = 1:Nr
    if isempty(current{i,1})
        [~,indCon] = sort(Cosine(:,i),'descend');
        k = 1;
        if length(indCon) > Nr
            while Indflag(indCon(k),1) == 1
                k=k+1;
            end
            Best(i,1) = indCon(k);
            Indflag(indCon(k),1) = 1;
        else
            Best(i,1) = indCon(1);
        end
    end
end
end
