function [best_x,fval]=BFO1cluster(fun,bacterialnum,iter_max,ped,LB,UB)%改进细菌觅食算法,聚类
[n,c]=size(LB);     %n类，c维
for i=1:bacterialnum
    bacterial(i).center=LB+(UB-LB).*rand(n,c);  %每类的聚类中心
    fitness(i,1)=fun(bacterial(i).center);
end
[a,b]=sort(fitness,'ascend');  %适应度函数及排序
bacterial=bacterial(b);
fval=fitness(1,1);
best_x=bacterial(1).center;  %全局最优解，求最小值
for iter=1:iter_max
    F=0.3*(exp(iter_max/(iter_max+iter))-1);   %复制
    b1=bacterialnum*0.5+ceil(bacterialnum*0.5/2);
    bacterial(bacterialnum*0.5+1:b1)=mutation_DE1(bacterial(bacterialnum*0.5+1:b1),F,[],LB,UB,1);
    for i=b1+1:bacterialnum
        bacterial(i).center=best_x;  
    end
    for i=1:bacterialnum
        if ped>rand
           bacterial(i).center=0.8.*best_x+rand(n,c);
           bacterial(i).center=boundtest(bacterial(i).center,LB,UB);
        end
    end
    for i=1:bacterialnum
       fitness(i,1)=fun(bacterial(i).center);
    end
    for i=1:ceil(bacterialnum*0.4)   %前40%一般
        %b2=10^((iter_max-iter)/iter_max);
        %step=(mysqrt(fitness(i,1),3)+0.6*b2)/(mysqrt(fitness(i,1),3)+30);
        step=(UB-LB)./(2*iter);
        temp_f=fitness(i,1);
        delta=unifrnd(-1,1,c,n);
        x4=bacterial(i).center+step.*delta'./norm(delta);
        x4=boundtest(x4,LB,UB);
        y4=fun(x4);
        if y4<temp_f
            bacterial(i).center=x4;
            fitness(i,1)=y4;
        end
    end
    bacterial(ceil(bacterialnum*0.4)+1:bacterialnum)=mutation_DE1(bacterial(ceil(bacterialnum*0.4)+1:bacterialnum),F,[],LB,UB,1);   %后60%差分
    for i=1:bacterialnum
        fitness(i,1)=fun(bacterial(i).center);
    end
    [a,b]=sort(fitness,'ascend');
    bacterial=bacterial(b);
    if a(1)<fval
       best_x=bacterial(1).center;
       fval=a(1);
    end   
end










