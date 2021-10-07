function [best_x,fval]=GWO(fun,popsize,iter_max,LB,UB)    %ª“¿«À„∑®
NC=size(LB,1);
for i=1:popsize
    wolf(i,:)=LB'+(UB-LB)'.*rand(1,NC);
    fitness(i,1)=fun(wolf(i,:));
end
[a,b]=sort(fitness,'ascend');
wolf=wolf(b,:);
alpha=wolf(1,:);
beta=wolf(2,:);
delta=wolf(3,:);
fval=a(1);
best_x=wolf(b(1),:);
for iter=1:iter_max
    a=2*cos((iter/iter_max)*pi/2);
    A=2*a.*rand(3,NC)-a;
    c=2.*rand(3,NC);
    new=alpha;
    for i=1:NC
        if rand<=1/NC
            new(i)=LB(i,1)+(UB(i,1)-LB(i,1))*rand;
        end
    end
    new=boundtest(new,LB,UB,2);
    newf=fun(new);
    if newf<fval
        fval=newf;
        best_x=new;
        wolf(1,:)=new;
        a(1)=newf;
        alpha=new;
    end  
    for i=1:popsize
        da=norm(c(1,:).*alpha-wolf(i,:));
        db=norm(c(2,:).*beta-wolf(i,:));
        dd=norm(c(3,:).*delta-wolf(i,:));
        x1=alpha-A(1,:).*da;
        x2=beta-A(2,:).*db;
        x3=delta-A(3,:)*dd;
        new=wolf(i,:);
        for j=1:NC
            if i==1
               k=ceil(NC*rand);
               if j==k
                   num=ceil(popsize*rand(1,2));
                   if num(1)==num(2)||num(1)==1||num(2)==1
                       num=ceil(popsize*rand(1,2));
                   end
                   new(j)=alpha(j)+2*a*rand*(wolf(num(1),j)-wolf(num(2),j));
               else
                   new(j)=alpha(k);
               end
            elseif i==2
                 if rand>0.67
                    new(j)=x1(j);
                 else
                     k=ceil(NC*rand);
                     new(j)=(alpha(k)+beta(k))/2;
                 end
            elseif i==3
                if rand>0.33
                    new(j)=(x1(j)+x2(j))/2;
                else
                    k=ceil(NC*rand); 
                    new(j)=(alpha(k)+beta(k)+delta(k))/3;
                end
            else
                new(j)=(x1(j)+x2(j)+x3(j))/3;
            end
        end
        new=boundtest(new,LB,UB,2);
        newf=fun(new);
        if newf<=fitness(i,1)
            wolf(i,:)=new;
            fitness(i,1)=newf;
        end
    end
    [a,b]=sort(fitness,'ascend');
    wolf=wolf(b,:);
    alpha=wolf(1,:);
    beta=wolf(2,:);
    delta=wolf(3,:);
    if a(1)<fval
       fval=a(1);
       best_x=wolf(1,:);
    end
end
    
   
    
    
        
    
