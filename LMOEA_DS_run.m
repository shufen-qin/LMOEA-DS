%% The code is written by Shufen Qin (This is a main file.)
function LMOEA_DS_run
clear
clc
%% Parameter settings
FEmax = 80000;
N = 153;
m = 3;
Ns = 30; % The number of solutions to be randomly generated along each guiding direction
problemLen =9;
functionvalue = cell(problemLen,1);
Dset = [500 1000 2000 5000];
for problemNum = 1:9
    problemNum
    problem = ['LSMOP',num2str(problemNum)];
    for i = 1:4
        D = Dset(i)
        PF = LSMOP('true',problem,m,10000,D);
        
        for  times= 1:20
            times
            %% Initialization
            [PopX,Boundary] = LSMOP('init',problem,m,N,D);
            PopX = unique(PopX,'rows');
            
            % Evaluation
            Obj = LSMOP('value',problem,m,PopX,D);
            FE = size(Obj,1);
            Pop = [PopX,Obj];
            
            % Reference vectors generation
            [W,N] = UniformPoint(N,m);
            W(W == 0) = 0.000001;
            
            %% Optimization
            while FE < FEmax
                [Arc,FE] = DirectedSampling(Pop,Ns,W,FE,Boundary,problem);
                [Pop,FE] = ComplementaryEnvironmentalSelection(Pop,Arc,W,FE,problem,Boundary);
            end
            
            %% save data
            functionvalue{problemNum,times} = Pop;
            IGDvalue(times,problemNum) = IGD(Pop(:,D+1:D+m),PF);
%             HVvalue(times,problemNum) = HV(Pop(:,D+1:D+m),PF);
            IGDmedian(problemNum) = median(IGDvalue(:,problemNum));
%             HVmedian(problemNum) = median(HVvalue(:,problemNum));
            IGDmean(problemNum) = mean(IGDvalue(:,problemNum));
%             HVmean(problemNum) = mean(HVvalue(:,problemNum));
            save(['Data\LMOEA-DS_',num2str(m),'m',num2str(D),'D.mat'],'IGDvalue','functionvalue','IGDmedian','IGDmean');
        end
    end
end
end
