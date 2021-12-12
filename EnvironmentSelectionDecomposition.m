function [Population] = EnvironmentSelectionDecomposition(Population,associate,Cosinemax,m)
% The environmental selection of the decomposition-based method.

% [Nr,m] = size(W);
[np,n] = size(Population);
f = Population(:,n-m+1:n);

%% Translate the population
f = (f-repmat(min(f),np,1))./(repmat(max(f),np,1)-repmat(min(f),np,1));

%% Associate each solution to a reference vector
% Cosinetemp = pdist2(f,W,'cosine');
% Cosine = 1-Cosinetemp;
% [Cosinemax,associate] = max(Cosine,[],2);
% Next = zeros(Nr,1);

%% Select one solution for each reference vector
list = unique(associate)';
Next = zeros(length(list),1);
t = 1;
for i = list
    current = find(associate == i);
    dist = pdist2(f(current,:),zeros(1,m),'Euclidean');
    Fan = Cosinemax(current)./dist;
    [~,best] = max(Fan);
    Next(t)  = current(best);
    t = t+1;
end
% Population for next generation
Population = Population(Next,:);

end
