function [best_x,fval]=IWO(fun,popsize,iter_max,LB,UB)   %野草算法
NC=size(LB,1);
popsize_max=ceil(2.0*popsize);
seed_min=1;
seed_max=10;
sigma_int=max(max(abs(LB)),max(abs(UB)))/3;
sigma_final=0.1;
for i=1:popsize
    weed(i,:)=LB'+(UB-LB)'.*rand(1,NC);
    fitness(i,1)=fun(weed(i,:));
end
[a,b]=sort(fitness,'ascend');
best_x=weed(b(1),:);
fval=a(1);
best=a(1);worst=a(end);
pop=weed;y=fitness;
for iter=1:iter_max
    sigma=(iter_max/iter-1)^3*(sigma_int-sigma_final)+sigma_final;
    for i=1:popsize
        seed_num(i)=floor((fitness(i)-worst)*(seed_max-seed_min)/(best-worst))+seed_min;
        for j=1:seed_num(i)
            new=boundtest(weed(i,:)+unifrnd(-sigma,sigma),LB,UB);
            pop=[pop;new];
            y=[y;fun(new)];
        end
    end
    num=size(pop,1);  %总的种群数
    if num>popsize_max
       [a1,b1]=sort(y,'ascend');  %排序
       pop=pop(b1,:);
       y=y(b1,1);
       weed=redu(pop,popsize_max+1:num,'r');
       fitness=redu(y,popsize_max+1:num,'r');
       a1=redu(a1,popsize_max+1:num,'r');
       b1=redu(b1,popsize_max+1:num,'r');
       if a1(1)<fval
           fval=y(1);
           best_x=pop(1,:);
       end
       best=y(1,1);worst=y(end,1);
       clear pop y
       popsize=popsize_max;
       pop=weed;y=fitness;
    else
       [a2,b2]=sort(y,'ascend');
       best=a2(1);worst=a2(end);
       if a2(1)<fval
          fval=a2(1);
          best_x=pop(b2(1),:);
       end
       clear a2 b2 
    end
end
   
    
    
        
        
