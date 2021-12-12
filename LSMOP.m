function [Output,Boundary] = LSMOP(Operation,Problem,M,Input,D)
% <Problem> <LSMOP>
% Large-scale benchmark MOP
% nk --- 5 --- Number of subcomponents in each variable group

%------------------------------- Reference --------------------------------
% R. Cheng, Y. Jin, and M. Olhofer, Test Problems for large-scale
% multiobjective and many-objective optimization, IEEE Transactions on
% Cybernetics, 2017, 47(12): 4108-4121.
%--------------------------------------------------------------------------


nk = 5; % Number of subcomponents in each variable group

lower    = zeros(1,D);
upper    = [ones(1,M-1),10.*ones(1,D-M+1)];
Boundary = [upper; lower ];

% Calculate the number of variables in each subcomponent
c = 3.8*0.1*(1-0.1);
for i = 1 : M-1
    c = [c,3.8.*c(end).*(1-c(end))];
end
sublen = floor(c./sum(c).*D/nk);
len    = [0,cumsum(sublen*nk)];
switch Operation
    case 'init'
        PopDec = rand(Input,D);
        PopDec = PopDec.*repmat(Boundary(2,:),Input,1)+(1-PopDec).*repmat(Boundary(1,:),Input,1);
        Output   = PopDec;
        
    case 'value'
        PopDec = Input;
        [N,D] = size(PopDec);
        G = zeros(N,M);
        switch Problem
            case 'LSMOP1'
                PopDec(:,M:D) = (1+repmat((M:D)./D,N,1)).*PopDec(:,M:D) - repmat(PopDec(:,1)*10,1,D-M+1);
                for i = 1 : 2 : M
                    for j = 1 : nk
                        G(:,i) = G(:,i) + Sphere(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                    end
                end
                for i = 2 : 2 : M
                    for j = 1 : nk
                        G(:,i) = G(:,i) + Sphere(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                    end
                end
                G      = G./repmat(sublen,N,1)./nk;
                PopObj = (1+G).*fliplr(cumprod([ones(N,1),PopDec(:,1:M-1)],2)).*[ones(N,1),1-PopDec(:,M-1:-1:1)];
                
            case 'LSMOP2'
                PopDec(:,M:D) = (1+repmat((M:D)./D,N,1)).*PopDec(:,M:D) - repmat(PopDec(:,1)*10,1,D-M+1);
                for i = 1 : 2 : M
                    for j = 1 : nk
                        G(:,i) = G(:,i) + Griewank(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                    end
                end
                for i = 2 : 2 : M
                    for j = 1 : nk
                        G(:,i) = G(:,i) + Schwefel(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                    end
                end
                G      = G./repmat(sublen,N,1)./nk;
                PopObj = (1+G).*fliplr(cumprod([ones(N,1),PopDec(:,1:M-1)],2)).*[ones(N,1),1-PopDec(:,M-1:-1:1)];
            case 'LSMOP3'
                PopDec(:,M:D) = (1+repmat((M:D)./D,N,1)).*PopDec(:,M:D) - repmat(PopDec(:,1)*10,1,D-M+1);
                for i = 1 : 2 : M
                    for j = 1 : nk
                        G(:,i) = G(:,i) + Rastrigin(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                    end
                end
                for i = 2 : 2 : M
                    for j = 1 : nk
                        G(:,i) = G(:,i) + Rosenbrock(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                    end
                end
                G      = G./repmat(sublen,N,1)./nk;
                PopObj = (1+G).*fliplr(cumprod([ones(N,1),PopDec(:,1:M-1)],2)).*[ones(N,1),1-PopDec(:,M-1:-1:1)];
                
            case 'LSMOP4'
                PopDec(:,M:D) = (1+repmat((M:D)./D,N,1)).*PopDec(:,M:D) - repmat(PopDec(:,1)*10,1,D-M+1);
                for i = 1 : 2 : M
                    for j = 1 : nk
                        G(:,i) = G(:,i) + Ackley(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                    end
                end
                for i = 2 : 2 : M
                    for j = 1 : nk
                        G(:,i) = G(:,i) + Griewank(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                    end
                end
                G      = G./repmat(sublen,N,1)./nk;
                PopObj = (1+G).*fliplr(cumprod([ones(N,1),PopDec(:,1:M-1)],2)).*[ones(N,1),1-PopDec(:,M-1:-1:1)];
                
            case 'LSMOP5'
                PopDec(:,M:D) = (1+repmat(cos((M:D)./D*pi/2),N,1)).*PopDec(:,M:D) - repmat(PopDec(:,1)*10,1,D-M+1);
                for i = 1 : 2 : M
                    for j = 1 : nk
                        G(:,i) = G(:,i) + Sphere(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                    end
                end
                for i = 2 : 2 : M
                    for j = 1 : nk
                        G(:,i) = G(:,i) + Sphere(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                    end
                end
                G      = G./repmat(sublen,N,1)./nk;
                PopObj = (1+G+[G(:,2:end),zeros(N,1)]).*fliplr(cumprod([ones(N,1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(N,1),sin(PopDec(:,M-1:-1:1)*pi/2)];
                
            case 'LSMOP6'
                PopDec(:,M:D) = (1+repmat(cos((M:D)./D*pi/2),N,1)).*PopDec(:,M:D) - repmat(PopDec(:,1)*10,1,D-M+1);
                for i = 1 : 2 : M
                    for j = 1 : nk
                        G(:,i) = G(:,i) + Rosenbrock(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                    end
                end
                for i = 2 : 2 : M
                    for j = 1 : nk
                        G(:,i) = G(:,i) + Schwefel(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                    end
                end
                G      = G./repmat(sublen,N,1)./nk;
                PopObj = (1+G+[G(:,2:end),zeros(N,1)]).*fliplr(cumprod([ones(N,1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(N,1),sin(PopDec(:,M-1:-1:1)*pi/2)];
                
            case 'LSMOP7'
                PopDec(:,M:D) = (1+repmat(cos((M:D)./D*pi/2),N,1)).*PopDec(:,M:D) - repmat(PopDec(:,1)*10,1,D-M+1);
                for i = 1 : 2 : M
                    for j = 1 : nk
                        G(:,i) = G(:,i) + Ackley(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                    end
                end
                for i = 2 : 2 : M
                    for j = 1 : nk
                        G(:,i) = G(:,i) + Rosenbrock(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                    end
                end
                G      = G./repmat(sublen,N,1)./nk;
                PopObj = (1+G+[G(:,2:end),zeros(N,1)]).*fliplr(cumprod([ones(N,1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(N,1),sin(PopDec(:,M-1:-1:1)*pi/2)];
                
            case 'LSMOP8'
                PopDec(:,M:D) = (1+repmat(cos((M:D)./D*pi/2),N,1)).*PopDec(:,M:D) - repmat(PopDec(:,1)*10,1,D-M+1);
                for i = 1 : 2 : M
                    for j = 1 : nk
                        G(:,i) = G(:,i) + Griewank(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                    end
                end
                for i = 2 : 2 : M
                    for j = 1 : nk
                        G(:,i) = G(:,i) + Sphere(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                    end
                end
                G      = G./repmat(sublen,N,1)./nk;
                PopObj = (1+G+[G(:,2:end),zeros(N,1)]).*fliplr(cumprod([ones(N,1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(N,1),sin(PopDec(:,M-1:-1:1)*pi/2)];
                
            case 'LSMOP9'
                PopDec(:,M:D) = (1+repmat(cos((M:D)./D*pi/2),N,1)).*PopDec(:,M:D) - repmat(PopDec(:,1)*10,1,D-M+1);
                for i = 1 : 2 : M
                    for j = 1 : nk
                        G(:,i) = G(:,i) + Sphere(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                    end
                end
                for i = 2 : 2 : M
                    for j = 1 : nk
                        G(:,i) = G(:,i) + Ackley(PopDec(:,len(i)+M-1+(j-1)*sublen(i)+1:len(i)+M-1+j*sublen(i)));
                    end
                end
                G = 1 + sum(G./repmat(sublen,N,1)./nk,2);
                PopObj(:,1:M-1) = PopDec(:,1:M-1);
                PopObj(:,M)     = (1+G).*(M-sum(PopObj(:,1:M-1)./(1+repmat(G,1,M-1)).*(1+sin(3*pi.*PopObj(:,1:M-1))),2));
                
        end
        Output = PopObj;
    case 'true'
        switch Problem
            case {'LSMOP1','LSMOP2','LSMOP3','LSMOP4'}
                P = UniformPoint(Input,M);
            case {'LSMOP5','LSMOP6','LSMOP7','LSMOP8'}
                P = UniformPoint(Input,M);
                P = P./repmat(sqrt(sum(P.^2,2)),1,M);
            case 'LSMOP9'
                interval     = [0,0.251412,0.631627,0.859401];
                median       = (interval(2)-interval(1))/(interval(4)-interval(3)+interval(2)-interval(1));
                X            = ReplicatePoint(Input,M-1);
                X(X<=median) = X(X<=median)*(interval(2)-interval(1))/median+interval(1);
                X(X>median)  = (X(X>median)-median)*(interval(4)-interval(3))/(1-median)+interval(3);
                P            = [X,2*(M-sum(X/2.*(1+sin(3*pi.*X)),2))];
                
        end
        Output = P;
end

end

function f = Sphere(x)
f = sum(x.^2,2);
end
function f = Griewank(x)
f = sum(x.^2,2)./4000 - prod(cos(x./repmat(sqrt(1:size(x,2)),size(x,1),1)),2) + 1;
end

function f = Schwefel(x)
f = max(abs(x),[],2);
end
function f = Rastrigin(x)
f = sum(x.^2-10.*cos(2.*pi.*x)+10,2);
end

function f = Rosenbrock(x)
f = sum(100.*(x(:,1:size(x,2)-1).^2-x(:,2:size(x,2))).^2+(x(:,1:size(x,2)-1)-1).^2,2);
end
function f = Ackley(x)
f = 20-20.*exp(-0.2.*sqrt(sum(x.^2,2)./size(x,2)))-exp(sum(cos(2.*pi.*x),2)./size(x,2))+exp(1);
end
function W = ReplicatePoint(SampleNum,M)
if M > 1
    SampleNum = (ceil(SampleNum^(1/M)))^M;
    Gap       = 0:1/(SampleNum^(1/M)-1):1;
    eval(sprintf('[%s]=ndgrid(Gap);',sprintf('c%d,',1:M)))
    eval(sprintf('W=[%s];',sprintf('c%d(:),',1:M)))
else
    W = (0:1/(SampleNum-1):1)';
end
end
