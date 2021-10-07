function [zbest_x,zbest_f]=psoBFO(fun,bacterialnum,iter_max,ped,LB,UB)%粒子群细菌觅食算法,求极大
NC=size(LB,1);
for i=1:bacterialnum
    x1=LB'+(UB-LB)'.*rand(1,NC);
    for j=1:NC
        if x1(j)>(LB(j,1)+UB(j,1))/2
            x2(j)=unifrnd((LB(j,1)+UB(j,1))/2,x1(j));
        else
            x2(j)=unifrnd(x1(j),(LB(j,1)+UB(j,1))/2);  
        end
    end
    if fun(x1)>fun(x2)     
       bacterial(i,:)=x1;
       fitness(i,1)=fun(x1);
    else
       bacterial(i,:)=x2;
       fitness(i,1)=fun(x2);
    end
    v(i,:)=rand(1,NC);
end
vmax=1;vmin=-1;
[a,b]=sort(fitness,'descend');
bacterial=bacterial(b,:);
fitness=fitness(b,1);
zbest_x=bacterial(1,:);
zbest_f=fitness(1,1);
pbest_x=bacterial;
pbest_f=fitness;
for iter=1:iter_max
    F=0.3*(exp(iter_max/(iter_max+iter))-1);   %复制
    b1=bacterialnum*0.5+ceil(bacterialnum*0.5/2);
    bacterial(bacterialnum*0.5+1:b1,:)=mutation_DE(bacterial(bacterialnum*0.5+1:b1,:),F,[],[LB UB],1);
    bacterial(bacterialnum*0.5+1:b1,:)=boundtest(bacterial(bacterialnum*0.5+1:b1,:),LB,UB);
    for i=b1+1:bacterialnum
        bacterial(i,:)=zbest_x;  
    end
    for i=1:bacterialnum
        if ped>rand
           bacterial(i,:)=0.8.*zbest_x+rand(1,NC);
           bacterial(i,:)=boundtest(bacterial(i,:),LB,UB);
        end
    end
    for i=1:bacterialnum
       fitness(i,1)=fun(bacterial(i,:));
       if fitness(i,1)>pbest_f(i,1)
           pbest_f(i,1)=fitness(i,1);
           pbest_x(i,:)=bacterial(i,:);
       end
    end
    [a,b]=sort(fitness,'descend');
    if a(1)>zbest_f
        zbest_f=a(1);
        zbest_z=bacterial(b(1),:);
    end
    bacterial=bacterial(b,:);
    w=0.9-iter*0.6/iter_max;
    for i=1:ceil(bacterialnum*0.4)   %前40%一般
        v(i,:)=w.*v(i,:)+2.05*rand.*(pbest_x(i,:)-bacterial(i,:))+2.05*rand.*(zbest_x-bacterial(i,:));
        v(i,find(v(i,:)>vmax))=vmax;
        v(i,find(v(i,:)<vmin))=vmin;
        bacterial(i,:)=bacterial(i,:)+v(i,:);
        bacterial(i,:)=boundtest(bacterial(i,:),LB,UB);
    end
    bacterial(ceil(bacterialnum*0.4)+1:bacterialnum,:)=mutation_DE(bacterial(ceil(bacterialnum*0.4)+1:bacterialnum,:),F,[],[LB UB],1);   %后60%差分
    bacterial(ceil(bacterialnum*0.4)+1:bacterialnum,:)=boundtest(bacterial(ceil(bacterialnum*0.4)+1:bacterialnum,:),LB,UB);
    for i=1:bacterialnum
        if rand>0.6
           LBi=min(bacterial,[],1);
           UBi=max(bacterial,[],1);
           for j=1:NC
              if LBi(j)<LB(j,1)
                 LBi(j)=LB(j);
              end
              if UBi(j)>UB(j,1)
                 UBi(j)=UB(j);
              end
              if bacterial(i,j)>(LBi(j)+UBi(j))/2
                 x2(j)=unifrnd((LBi(j)+UBi(j))/2,bacterial(i,j));
              else
                 x2(j)=unifrnd(bacterial(i,j),(LBi(j)+UBi(j))/2);  
              end
           end
           if fun(x2)>fun(bacterial(i,:))     
              bacterial(i,:)=x2;
           end
           fitness(i,1)=fun(bacterial(i,:));
           if fitness(i,1)>pbest_f(i,1)
               pbest_x(i,:)=bacterial(i,:);
               pbest_f(i,1)=fitness(i,1);
           end
        end
    end
    [a,b]=sort(fitness,'descend');
    bacterial=bacterial(b,:);
    if a(1)>zbest_f
       zbest_x=bacterial(1,:);
       zbest_f=a(1);
    end   
end










