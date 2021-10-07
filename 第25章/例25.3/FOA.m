function [best_x,fval]=FOA(fun,popsize,iter_max,varargin)    %果蝇算法
type=varargin{end};
if type==1
    NC=varargin{1};
    X_axis=1e-6*rand;
    Y_axis=1e-6*rand;
    X=zeros(popsize,NC);
    Y=zeros(popsize,NC);
    D=zeros(popsize,NC);
    x=zeros(1,NC);
    for i=1:popsize
        X(i,:)=X_axis+2*rand(1,NC)-1;
        Y(i,:)=Y_axis+2*rand(1,NC)-1;
        for j=1:NC
            D(i,j)=sqrt(X(i,j)^2+Y(i,j)^2);
            %D(i,j)=abs(X(i,j))+abs(Y(i,j));
            S(i,j)=1/D(i,j);
            x(j)=S(i,j);
        end  
        smell(i,1)=fun(x);
    end
    [bestsmell,bestindex]=min(smell);
    X_axis=X(bestindex,:);
    Y_axis=Y(bestindex,:);
    best_x=S(bestindex,:);
    smellbest=bestsmell;
    for iter=1:iter_max
       for i=1:popsize
           X(i,:)=X_axis+2*rand(1,NC)-1;
           Y(i,:)=Y_axis+2*rand(1,NC)-1;
           for j=1:NC
               D(i,j)=sqrt(X(i,j)^2+Y(i,j)^2);
              % D(i,j)=abs(X(i,j))+abs(Y(i,j));
               S(i,j)=1/D(i,j);
               x(j)=S(i,j);
           end
           smell(i,1)=fun(x);
       end
       [bestsmell,bestindex]=min(smell);
       if bestsmell<smellbest
           X_axis=X(bestindex,:);
           Y_axis=Y(bestindex,:);
           smellbest=bestsmell;
           best_x=S(bestindex,:);
           fval=bestsmell;
       end
       yy(iter)=smellbest;
       xbest(iter,:)=X_axis;
       ybest(iter,:)=Y_axis;
    end
    figure(1)
    plot(yy);
    title('优化过程');
    xlabel('迭代数');ylabel('气味');
    figure(2)
    plot(xbest,ybest,'b.');
    title('果蝇飞行路径');
    xlabel('X-axis');ylabel('Y-axis');
elseif type==2   %处理负数
    LB=varargin{1};
    UB=varargin{2};
    NC=size(LB,1);
    R=1;
    X_axis=LB'+(UB-LB)'.*rand(1,NC);
    for i=1:popsize
        X(i,:)=X_axis+R.*unifrnd(-1,1,1,NC);
        X(i,:)=boundtest(X(i,:),LB,UB);
        smell(i)=fun(X(i,:));
    end
    [bestsmell,bestindex]=min(smell);
    X_axis=X(bestindex,:);
    smellbest=bestsmell;
    for iter=1:iter_max
        R=R*0.85^iter;
        for i=1:popsize
           X(i,:)=X_axis+R.*unifrnd(-1,1,1,NC);
           X(i,:)=boundtest(X(i,:),LB,UB);
           smell(i)=fun(X(i,:));
        end
        [bestsmell,bestindex]=min(smell);
        if bestsmell<smellbest
           X_axis=X(bestindex,:);
           smellbest=bestsmell;
           best_x=X(bestindex,:);
           fval=bestsmell;
        end
    end
end
    
    
    
    



     
