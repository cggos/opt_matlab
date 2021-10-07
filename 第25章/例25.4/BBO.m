function [best_x,fval]=BBO(fun,popsize,iter_max,LB,UB)   %生物地理算法
NC=size(LB,1);
I=1;E=1;
m_max=0.01;
k=ceil(0.1*popsize);    %精英个体
for i=1:popsize
    habitat(i,:)=LB'+(UB-LB)'.*rand(1,NC);
    fitness(i,1)=fun(habitat(i,:));
end
[a,b]=sort(fitness,'ascend');
best_x=habitat(b(1),:);
fval=a(1);
fmax=a(1);
habitat=habitat(b,:);
favg=mean(fitness);
fitness=fitness(b,1);
flag=0;
for iter=1:iter_max
    for i=1:popsize
        lamda(i,1)=I*(cos(pi*fitness(i,1)/fmax)+1)/2;
     %   lamda(i,1)=I*(1-(popsize-i)/popsize);
        mu(i,1)=E*(-cos(pi*fitness(i,1)/fmax)+1)/2;
     %  mu(i,1)=E*(popsize-i)/popsize;
    end
    for i=1:popsize
        if rand<1
           for j=1:NC
               if rand<lamda(i,1)  %需迁移
                  for k=1:popsize
                      if k~=i
                          if rand<mu(k,1)    %迁出的栖息地
                             habitat(i,j)=rand*habitat(i,j)+(1-rand)*habitat(k,j);
                             break  
                          end
                      end
                  end  
               end
           end
        end
        if i>k+1
            mi=m_max*abs(fitness(i,1)-favg)/favg;
            for j=1:NC
                  if rand<mi
                      habitat(i,:)=LB'+(UB-LB)'.*rand(1,NC);
                  end
             end
         end
         habitat(i,:)=boundtest(habitat(i,:),LB,UB);
         fitness(i,1)=fun(habitat(i,:));
    end
    [a,b]=sort(fitness,'ascend');
    favg=mean(fitness);
    habitat=habitat(b,:);
    fitness=fitness(b,1);
    fmax=a(1);
    if a(1)<fval
        fval=a(1);
        best_x=habitat(1,:);
        flag=0;
    else
        flag=flag+1;
    end
    if flag==3
        new=habitat(1,:)*(1+rand);
        if fun(new)<fval
           habitat(1,:)=new;
           fitness(1,1)=fun(new);
           fmax=fun(new);
           flag=0;
        end
    end 
end
    
        
        
                
