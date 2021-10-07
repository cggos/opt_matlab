function [best_x,best_f]=BFO1(fun,bacterialnum,iter_max,ped,LB,UB)%改进细菌觅食算法,求极大
NC=size(LB,1);
for i=1:bacterialnum
    for j=1:NC
        x0=rand;
        while x0==0.25||x0==0.5||x0==0.75
           x0=rand;
        end
        for k=1:300
           x0=4.*x0.*(1-x0);
        end
        x1(i,j)=LB(j,1)+(UB(j,1)-LB(j,1))*x0;
    end
    y1(i,1)=fun(x1(i,:));
end
for i=1:bacterialnum
    x2(i,:)=LB'+UB'-x1(i,:);
    y2(i,1)=fun(x2(i,:));
end
x3=[x1;x2];
y3=[y1;y2];
[a,b]=sort(y3,'descend');
x3=x3(b,:);
y3=y3(b,1);
for i=1:bacterialnum
    bacterial(i,:)=x3(i,:);
    fitness(i,1)=y3(i,1);
end
best_x=bacterial(1,:);
best_f=fitness(1,1);
for iter=1:iter_max
    F=0.3*(exp(iter_max/(iter_max+iter))-1);   %复制
    b1=bacterialnum*0.5+ceil(bacterialnum*0.5/2);
    bacterial(bacterialnum*0.5+1:b1,:)=mutation_DE(bacterial(bacterialnum*0.5+1:b1,:),F,[],[LB UB],1);
    bacterial(bacterialnum*0.5+1:b1,:)=boundtest(bacterial(bacterialnum*0.5+1:b1,:),LB,UB);
    for i=b1+1:bacterialnum
        bacterial(i,:)=best_x;  
    end
    for i=1:bacterialnum
        if ped>rand
           bacterial(i,:)=0.8.*best_x+rand(1,NC);
           bacterial(i,:)=boundtest(bacterial(i,:),LB,UB);
        end
    end
    for i=1:bacterialnum
       fitness(i,1)=fun(bacterial(i,:));
    end
    for i=1:ceil(bacterialnum*0.4)   %前40%一般
        %b2=10^((iter_max-iter)/iter_max);
        %step=(mysqrt(fitness(i,1),3)+0.6*b2)/(mysqrt(fitness(i,1),3)+30);
        step=(UB-LB)'./(2*iter);
        temp_f=fitness(i,1);
        delta=unifrnd(-1,1,1,NC);
        x4=bacterial(i,:)+step.*delta./norm(delta);
        x4=boundtest(x4,LB,UB);
        y4=fun(x4);
        if y4>temp_f
            bacterial(i,:)=x4;
            fitness(i,1)=y4;
        end
    end
    bacterial(ceil(bacterialnum*0.4)+1:bacterialnum,:)=mutation_DE(bacterial(ceil(bacterialnum*0.4)+1:bacterialnum,:),F,[],[LB UB],1);   %后60%差分
    bacterial(ceil(bacterialnum*0.4)+1:bacterialnum,:)=boundtest(bacterial(ceil(bacterialnum*0.4)+1:bacterialnum,:),LB,UB);
    for i=1:bacterialnum
        fitness(i,1)=fun(bacterial(i,:));
    end
    [a,b]=sort(fitness,'descend');
    bacterial=bacterial(b,:);
    if a(1)>best_f
       best_x=bacterial(1,:);
       best_f=a(1);
    end   
end










