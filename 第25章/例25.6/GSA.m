function [best_x,fval]=GSA(fun,popsize,iter_max,LB,UB)   %引力算法
G=100;alpha=20;beta=0.02;
NC=size(LB,1);
for i=1:popsize
    pop(i).x=LB'+(UB-LB)'.*rand(1,NC);
    pop(i).fitness=fun(pop(i).x);
    pop(i).v=rand(1,NC);
    y(i)=pop(i).fitness;
end
[a,b]=sort(y,'ascend');   %求极小
worst=a(end);xworst=pop(b(end)).x;
best=a(1);xbest=pop(b(1)).x;
cbest=pop(b(1));
pop=pop(b);
flag=0;
for iter=1:iter_max
    Gt=G*exp(-alpha*iter/iter_max);
    N=uint16((beta+(1-iter/iter_max)*(1-beta))*popsize);
    xworst=xworst.*(1+0.5*rand);
    xworst=boundtest(xworst,LB,UB);
    worst=fun(xworst);
    if flag==3
        xbest=xbest.*(1+0.5*rand);
    %    zbest=(xbest-LB')./(UB-LB)';
    %    for j=1:20
    %        zbest=4.0*zbest.*(1-zbest);
    %    end
    %    xbest=LB'+(UB-LB)'.*zbest;
        xbest=boundtest(xbest,LB,UB);
        best1=fun(xbest);
        if best1<best
            best=best1;
            cbest.x=xbest;
            cbest.fitness=best;
        end
    end
    M=0;
    for i=1:popsize
        pop(i).m=(pop(i).fitness-worst)/(best-worst);
        M=M+pop(i).m;
    end 
    for i=1:popsize
        pop(i).m=pop(i).m/M;
        k=0;
        for j=1:popsize
            if j~=i
                k=k+1;
                f(k,:)=Gt*pop(i).m*pop(j).m*(pop(j).x-pop(i).x)/(norm(pop(i).x-pop(j).x)+0.00001);
                f(k,:)=rand.*f(k,:);
            end
        end
        pop(i).F=sum(f(1:N-1,:));
        pop(i).a=pop(i).F./pop(i).m;
        pop(i).v=rand(1,NC).*pop(i).v+pop(i).a;
        pop(i).x=pop(i).x+pop(i).v;
        pop(i).x=boundtest(pop(i).x,LB,UB);
        pop(i).fitness=fun(pop(i).x);
        y(i)=pop(i).fitness;
    end
    [a,b]=sort(y,'ascend');   %求极小
    xworst=pop(b(end)).x;
    best=a(1);xbest=pop(b(1)).x;
    if a(1)<cbest.fitness
        cbest=pop(b(1));
        flag=0;
    else
        flag=flag+1;
    end
    pop=pop(b);
end
best_x=cbest.x;
fval=cbest.fitness;
               
                
                
                