function [bestx,fval]=SAPSO(fun,popsize,c1,c2,lamda,iter_max,LB,UB)
D=size(LB,1);
for i = 1:popsize
    x(i,:) =LB'+(UB-LB)'.*rand(1,D);  %随机初始化位置
    v(i,:) = randn(1,D);  %随机初始化速度
    p(i)=fun(x(i,:));
    y(i,:)=x(i,:);
end
[a,b]=min(p);
pg=x(b,:);             %Pg为全局最优
fval=a;
T=fun(pg)/log(5);   %第3步：定义初始温度
for iter = 1:iter_max
    groupFit = fun(pg);    
    for i = 1:popsize          %第4步：当前温度T下各个pi的适应值
        Tfit(i)=exp(-(p(i)-groupFit)/T);       
    end    
    SumTfit=sum(Tfit);    
    Tfit=Tfit/SumTfit;    
    pBest=rand();    
    for i=1:popsize          %第5步：采用轮盘赌策略确定全局最优的某个替代值        
        ComFit(i)=sum(Tfit(1:i));        
        if pBest<=ComFit(i)            
            pg_plus=x(i,:);            
            break;            
        end        
    end    
    C=c1+c2;    
    ksi=2/abs(2-C-sqrt(C^2-4*C));     %速度压缩因子   
    for i=1:popsize                               %更新各微粒的速度和位置
        v(i,:)=ksi*(v(i,:)+c1*rand*(y(i,:)-x(i,:))+c2*rand*(pg_plus-x(i,:)));
        x(i,:)=x(i,:)+v(i,:);
        if fun(x(i,:))<p(i)   %更新各微粒的pi值及群体的pg值
            p(i)=fun(x(i,:));
            y(i,:)=x(i,:);
        end
        if p(i)<fval
            pg=y(i,:);
            fval=p(i);
        end
    end
    T=T*lamda;       %降温   
end                 
bestx=pg;
