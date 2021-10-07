function [bestx,fval]=MA(fun,popsize,iter_max,Nc,a,b,LB,UB)   %猴群算法
NC=size(LB,1);
for i=1:popsize
    x1(i,:)=LB'+(UB-LB)'.*rand(1,NC);
    x2(i,:)=LB'+UB'-x1(i,:);    %相反点
    fitness1(i,1)=fun(x1(i,:));
    fitness2(i,1)=fun(x2(i,:));
end
pop=[x1;x2];
fitness=[fitness1;fitness2];
[a1,b1]=sort(fitness,'ascend');
monkey=pop(b1(1:popsize),:);
fitness=fitness(b1(1:popsize),1);
bestx=monkey(1,:);
fval=fitness(1,1);
flag=0;
N=ceil(popsize/2);
for iter=1:iter_max
    for i=1:N
        [monkey(i,:),fitness(i,1)]=climb(fun,monkey(i,:),Nc,iter,iter_max,a,b,LB,UB);   %爬
    end
    [a1,b1]=sort(fitness,'ascend');
    best=a1(1);worst=a1(end);
    temp=unifrnd(0.95,0.98);
    beta1=temp*log(iter)/log(iter_max);
    for i=N+1:popsize     %望跳
       beta2=(1-temp)*(worst-best+1)/(fitness(i,1)-best+1);
       nb=(beta1+beta2)*b;
       n=0;
       while n<=40     %望跳次数
            for j=1:NC
                LB1=monkey(i,j)-nb;
                if LB1<LB(j,1)
                    LB1=LB(j,1);
                end
                UB1=monkey(i,j)+nb;
                if UB1>UB(j,1)
                   UB1=UB(j,1);
                end
                new(j)=LB1+(UB1-LB1)*rand;
            end
            newF=fun(new);
            if newF<fitness(i,1)
               monkey(i,:)=new;
               fitness(i,1)=newF;
               break
            else
                n=n+1;
            end
       end
       [monkey(i,:),fitness(i,1)]=climb(fun,monkey(i,:),Nc,iter,iter_max,a,b,LB,UB);
    end
    [a1,b1]=sort(fitness,'ascend');
    worst=a1(end);
    best_x=monkey(b(1),:);
    [best_x,best]=PS(fun,best_x,1,1,1e-6,3000);   %猴王模式搜索
    for i=1:popsize     %翻
         if i~=b1(1)
            if rand>0.5
               v=1;
            else
               v=-1;
            end
            lamda=v*(log(fitness(i,1)-best+1)+1)/(log(worst-best+1)+1);
            for j=1:NC
                new(j)=monkey(i,j)+lamda*(best_x(j)-monkey(i,j));
            end
            new=boundtest(new,LB,UB);
            newF=fun(new);
            if newF<fitness(i,1)
               d=norm(new-best_x);
               if d>=0.05
                  monkey(i,:)=new;
                  fitness(i,1)=newF;
               else
                   monkey(i,:)=LB'+(UB-LB)'.*rand(1,NC);
                   fitness(i,1)=fun(monkey(i,:));
               end
            elseif newF<best
                best_x=new;
                best=newF;
                monkey(b1(1),:)=new;
                fitness(b1(1),1)=newF;
            end
         end
     end
     [a1,b1]=sort(fitness,'ascend');
     monkey=monkey(b1,:);
     fitness=fitness(b1,1);
     if a1(1)<fval
        flag=0;
        fval=a1(1);
        bestx=monkey(1,:);
     else
        flag=flag+1;
     end
     if flag==20   %高斯变异
        monkey(1,:)=bestx*(1+randn);
        monkey(1,:)=boundtest(monkey(1,:),LB,UB);
        fitness(1,1)=fun(monkey(1,:));
        if fitness(1,1)<fval
            fval=fitness(1,1);
            bestx=monkey(1,:);
        end
        flag=0;
    end
end
            
function  [pop,y]=climb(fun,pop,Nc,iter,iter_max,a,b,LB,UB)
NC=length(pop);
a0=b*sin((1-(iter-1)/iter_max)*pi/2);
for k=1:Nc
    na=a/10+a0*log(Nc/k)/log(Nc);
    for j=1:NC
        if rand>0.5
            delta(j)=na;
        else
            delta(j)=-na;
        end
        plusx(j)=pop(j)+delta(j);
        subx(j)=pop(j)-delta(j);
    end
    plusx=boundtest(plusx,LB,UB);
    subx=boundtest(subx,LB,UB);
    ch=fun(plusx)-fun(subx);
    for j=1:NC       
         F=ch/(2*delta(j));
         new(j)=pop(j)-na*sign(F);
    end
    if boundtest1(new,LB,UB)==1
       pop=new;
    end
end
y=fun(pop);

function [besty,f]=PS(fun,x0,d,alpha,esp,T)   %模式搜索
n=length(x0);
x(1,:)=x0;
y(1,:)=x(1,:);
h=eye(n);
k=1;
while k<T
  j=1;
  while 1
     a=y(j,:)+d*h(j,:);
     fa=fun(a);
     if fa<fun(y(j,:))
        y(j+1,:)=a;
     else
       a=y(j,:)-d*h(j,:); 
       fa=fun(a);
       if fa<fun(y(j,:))
          y(j+1,:)=a;
       else
          y(j+1,:)=y(j,:);  
       end
     end
     if j<n
         j=j+1;
     else
         break;
     end
  end
  if fun(y(n+1,:))<fun(x(k,:))
       x(k+1,:)=y(n+1,:);
       y(1,:)=x(k+1,:)+alpha*(x(k+1,:)-x(k,:));
       k=k+1;
  else
       if d<=esp
          besty=x(k,:);
          f=fun(besty);
          break;
       else
          y(1,:)=x(k,:);
          x(k+1,:)=x(k,:);
          d=d/2;
          k=k+1;
        end
  end
end
if exist('besty','var')==0
    besty=x0;
    f=fun(x0);
end
                