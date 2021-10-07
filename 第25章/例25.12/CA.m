function [bestx,fval]=CA(fun,popsize,iter_max,px,pm,LB,UB)  %文化算法
c=uint16(popsize*1.5);                   %从种群中随机取出的个体数
NC=size(LB,1);
for i=1:popsize
    pop(i).x=LB'+(UB-LB)'.*rand(1,NC);      %实数编码结构
    fitness(i,1)=fun(pop(i).x);
end
[fbest,bestnum]=min(fitness);
situation=pop(bestnum).x;            %初始化形势知识
for i=1:popsize
    for j=1:NC
        standard(i).l(j)=LB(j,1);              %规范知识
        standard(i).L(j)=inf;
        standard(i).u(j)=UB(j,1);              
        standard(i).U(j)=inf;
    end
end
for iter=1:iter_max
    acc=0.4*popsize*(1+1/iter);   %接受概率
    for i=1:popsize
        for j=1:NC
            if(pop(i).x(j)<situation(j))
                popson(i).x(j)=pop(i).x(j)+abs((standard(i).u(j)-standard(i).l(j))*rand);  %生成的子代个体
            elseif(pop(i).x(j)>situation(j))
                popson(i).x(j)=pop(i).x(j)-abs((standard(i).u(j)-standard(i).l(j))*rand);
            elseif(pop(i).x(j)==situation(j))
                popson(i).x(j)=pop(i).x(j)+(standard(i).u(j)-standard(i).l(j))*rand;
            end
        end
    end
  %  pop=[pop;popsize];
    for i=(popsize+1):(2*popsize)
        pop(i)=popson(i-popsize);     %将生成的子代个体放到父代个体后面   
    end
    winnum=[];            %每个个体的胜利次数
    for i=1:2*popsize     %选择个体
        cnum=randperm(2*popsize);
        winsingle=0;
        for j=1:c
            if(fun(pop(i).x)<fun(pop(cnum(j)).x)) %每个个体与随机的c个个体作比较
               winsingle=winsingle+1;          %记录胜利次数
            end
        end
        winnum(i)=winsingle; 
    end
    index=1;
    for i=1:2*popsize  
        for j=(i+1):2*popsize
            if(winnum(i)<winnum(j))    %将胜利次数从大到小排列
                uusee=winnum(i);
                winnum(i)=winnum(j);
                winnum(j)=uusee;           
                index=j;
            end       
        end
        pptv=pop(i);
        pop(i)=pop(index);   %将获胜次数多的前40个个体放入父代
        pop(index)=pptv;
    end
    for i=1:2*popsize                                                                                                                                 
        y(i,1)=fun(pop(i).x);       %计算新生成的父代个体的适应度值
    end
    [a,b]=sort(y,'ascend');
    pop=pop(b);
    y=y(b,1);
    pop=pop(1:popsize);
    fitness=y(1:popsize,1);
    pop=trun(pop);    %截尾
    if NC>5
       pop=back(pop,LB,UB);   %逆转
    end
    pop=cr(pop,px,LB,UB);    %交叉
    pop=m(pop,pm,LB,UB,iter,iter_max);    %变异
    for i=1:popsize
        fitness(i,1)=fun(pop(i).x);
    end
    [a,b]=sort(fitness,'ascend');
    fbest=a(1);bestnum=b(1);  %取出新生种群的最优个体
    if(fbest<fun(situation))
       if(rand<=acc)            %以概率acc择优接受形势知识
           situation=pop(bestnum).x;      
       end
    end
    for i=1:popsize
        for j=1:NC
           if pop(i).x(j)<standard(i).l(j)||fitness(i,1)<standard(i).L(j) 
              if(rand<=acc)            %以概率acc择优接受形势知识
                   standard(i).l(j)=pop(i).x(j);
                   standard(i).L(j)=fun(pop(i).x);
              end
           end
           if pop(i).x(j)>standard(i).u(j)||fitness(i,1)<standard(i).U(j) 
              if(rand<=acc)            %以概率acc择优接受形势知识
                  standard(i).u(j)=pop(i).x(j);
                  standard(i).U(j)=fun(pop(i).x);
              end
           end
        end
    end
end
bestx=pop(bestnum).x;
fval=fbest;

function pop=m(pop,pm,LB,UB,iterm,iterm_max)
popsize=size(pop,2);
genelen=length(pop(1).x);
for i=1:popsize
    for j=1:genelen
        if rand<pm 
           v1=pop(i).x(j)-UB(j);
           v2=LB(j)-pop(i).x(j);
           fg=rand*(1-iterm/iterm_max)^2;
           if rand>0.5
              pop(i).x(j)=pop(i).x(j)+v1*fg;
           else
              pop(i).x(j)=pop(i).x(j)+v2*fg;
           end
           pop(i).x(j)=boundtest(pop(i).x(j),LB(j),UB(j));
       end
   end
end 

function pop=trun(pop)   %截尾算子
popsize=size(pop,2);
num=uint16(0.4*popsize);
new=[pop(1:num) pop(1:num) pop];
pop=new(1:popsize);

function pop=back(pop,LB,UB)
popsize=size(pop,2);
NC=length(pop(1).x);
for i=1:popsize
    new=pop(i).x;
    point=sort(ceil(NC*rand(1,2)),'ascend');
    while point(1)==point(2)||point(2)-point(1)<2
        point=sort(ceil(NC*rand(1,2)),'ascend');
    end
    pop(i).x=[new(1:point(1)) fliplr(new(point(1)+1:point(2)-1)) new(point(2):end)];
    pop(i).x=boundtest(pop(i).x,LB,UB);
end


