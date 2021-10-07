function [best_x,fval]=APOA(fun,popsize,iter_max,LB,UB)    %人工植物算法
NC=size(LB,1);
alpha=-2.0011e-5;pmax=0.034;rd=0.57;
pm=0.5;theta=pi/10;w=pi/2;re=0.25;rate=0.8;
growthmin=max(UB,[],1)/50;growthmax=max(UB,[],1)/10;
angle=(2/pi).*ones(popsize,NC);
c=exp(-max(UB,[],1)/2)/sqrt(2*pi);
pxmax=0.80;
for i=1:popsize
    plant(i,:)=LB'+(UB-LB)'.*rand(1,NC);
    fitness(i,1)=fun(plant(i,:));
end
[a,b]=sort(fitness,'ascend');
best=a(1);worst=a(end);
best_x=plant(b(1),:);
fval=a(1);
for iter=1:iter_max
    if mod(iter,5)==0
        growth=growthmax;
    elseif mod(iter,5)~=0
        growth=growthmin;
    else
        growth=(growthmax-growthmin)*(iter/iter_max)^2+(growthmax-growthmin)*(iter/iter_max);
    end
    for i=1:popsize
        uf(i,1)=(worst-fitness(i,1))/(worst-best);
        %p(i)=alpha*uf(i,1)*pmax/(alpha*uf(i,1)+pmax)-rd;
        p(i)=alpha*uf(i,1)^2+pmax*uf(i,1)-rd;
    end
    [a,b]=sort(p,'descend');
    plant=plant(b,:);
    p=p(b);
    angle=angle(b,:);
    num=ceil(pm*popsize);
    bestx=plant(1,:);
    for i=1:popsize
        k=2*sqrt(abs(p(i))/(0.25*pi*re));
        for j=1:NC
           theta1=k*sqrt(abs(cos(theta+angle(i,j))-cos(w+angle(i,j))));
           angle(i,j)=angle(i,j)+theta1;
        end
        a1=1;
        for j=1:NC-1
            a1=a1*cos(angle(i,j));
        end
        for j=1:NC
            if j==1
                d(j)=a1;
            elseif j==NC
                d(j)=sin(angle(i,j-1));
            else
                d(j)=sin(angle(i,j))*a1;
            end
        end
        if i<=num
            plant(i,:)=plant(i,:)+(bestx-plant(i,:))*growth*rand;
        else
            plant(i,:)=plant(i,:)+growth*rand.*d; 
        end
        plant(i,:)=boundtest(plant(i,:),LB,UB);
        fitness(i,1)=fun(plant(i,:));
    end
    [a,b]=sort(fitness,'ascend');
    bestx=plant(b(1),:);
    favg=mean(fitness);
    if max(abs(fitness-favg))>=1
        f=1;
    else
        f=max(abs(fitness-favg));
    end
    sigma=sum(((fitness-favg)./f).^2)/popsize;
   % px=(pmax-pmin)*sigma/iter_max+(pmin-pmax)*2*sigma/iter_max+pmax;
    px=pxmax*(1-sigma/popsize);
    if rand<px
        for i=1:popsize
            plant(i,:)=plant(i,:)+(bestx-plant(i,:))*growth*rand+c*rand;
            plant(i,:)=boundtest(plant(i,:),LB,UB);
            fitness(i,1)=fun(plant(i,:));
        end
    end
    [a,b]=sort(fitness,'ascend');
    worst=a(end);
    bestx=plant(b(1),:);        
    if rand<rate
        bestx=bestx+(bestx-plant(b(end),:))*growth*rand;
    elseif rand<1/NC
        bestx=bestx+(UB-LB)'*rand;
    end
    bestx=boundtest(bestx,LB,UB);
    fitness(b(1),1)=fun(bestx);
    plant(b(1),:)=bestx;
    best=fitness(b(1),1);
    if fitness(b(1),:)<fval
        fval=fitness(b(1),1);
        best_x=bestx;
    end     
end
    
        
        
    
    
  