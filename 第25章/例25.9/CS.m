function [best_x,fval]=CS(fun,popsize,iter_max,LB,UB)  %布谷鸟算法
NC=size(LB,1);
for i=1:popsize
    bird(i,:)=LB'+(UB-LB)'.*rand(1,NC);
    fitness(i,1)=fun(bird(i,:));
end
[a,b]=min(fitness);
best_x=bird(b,:);
fval=a;
beta=1.5;
%a1=mfun('gamma',beta+1);  %levy步长
a1=gamma(beta+1);  %levy步长
%a2=mfun('gamma',(1+beta)/2);
a2=gamma((1+beta)/2);
a3=beta*2^((beta-1)/2);
sigmau=(a1*sin(pi*beta/2)/(a2*a3))^(1/beta);
for iter=1:iter_max
    pa=1-0.95*iter/iter_max;
    w=0.05+0.95*unifrnd(0,1)+0.01*randn;
    alpha=0.01+0.64*exp(-20*(iter/iter_max)^4);
    v=abs(normrnd(0,1));
    u=normrnd(0,sigmau^2);
    L=u/(v)^(1/beta);
    for i=1:popsize
        new=bird(i,:)+rand*(best_x-bird(i,:));
        new=boundtest(new,LB,UB);
        newF=fun(new);
        if newF<fitness(i,1)
            bird(i,:)=new;
            fitness(i,1)=newF;
        end
    end
    new1=bird;
    for i=1:popsize
        new1(i,:)=w*new1(i,:)+alpha*L.*(new1(i,:)-best_x);
        new1(i,:)=boundtest(new1(i,:),LB,UB);
        fitness1(i,1)=fun(new1(i,:));
    end
    pop=[bird;new1];
    y=[fitness;fitness1];
    [a4,b4]=sort(y,'ascend');
    pop=pop(b4,:);
    y=y(b4,1);
    bird=pop(1:popsize,:);
    fitness=y(1:popsize,1);
    for i=1:popsize
        if rand>pa
            new=bird(i,:);
            k=ceil(rand*NC);
            new(k)=LB(k,1)+(UB(k,1)-LB(k,1))*rand;
            f=fun(new);
            if f>fitness(i,1)
                new=bird(i,:);
            end
            new=new*(1+0.03*randn);
            new=boundtest(new,LB,UB);
            f=fun(new);
            if f<fitness(i,:)
                bird(i,:)=new;
                fitness(i,1)=f;
            end
        end
    end
    [a,b]=min(fitness);
    if a<fval
        best_x=bird(b,:);
        fval=a;
    end
end
            
            
            
    
        