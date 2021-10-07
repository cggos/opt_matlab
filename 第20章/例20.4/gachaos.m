function cbest=gachaos(fun,popsize,iterm_max,pm,px,LB,UB)   %遗传-混沌算法
if isempty(pm)
    pm=0.1;
end
if isempty(px)
    px=0.9;
end
if isempty(iterm_max)
    iterm_max=3000;
end
if isempty(popsize)
    popsize=30;
end
pop=in(popsize,LB,UB);
for i=1:popsize
    pop(i).fitness=fun(pop(i).x);
    y(i)=pop(i).fitness;
end
[a,b]=sort(y,'descend');
cbest=pop(b(1));
%some=ceil(popsize/3);  %随机选1/3个体作混沌搜索
for iterm=1:iterm_max
    pop=cal(pop);    %适应度值
    pop=sele(pop);     %选择
    pop=cr(pop,px,LB,UB); %交叉
    pop=mut(pop,pm,LB,UB,iterm,iterm_max); %变异
    %pop=pop(randperm(popsize));
    %somepop=pop(1:some);   
    %pop(1:some)=chaos_serch(fun,somepop,LB,UB);
    pop=chaos_serch(fun,pop,LB,UB);
    for i=1:popsize
        pop(i).fitness=fun(pop(i).x);
        y(i)=pop(i).fitness;
    end
    [a,b]=sort(y,'descend');
    if a(1)>cbest.fitness
        cbest=pop(b(1));
        pop(b(end))=pop(b(1));
    end
end

function pop=cal(pop)
popsize=size(pop,2);
for i=1:popsize
    pop(i).index=-1;
end
for i=1:popsize
    index=1;
    for j=1:popsize
        if pop(j).fitness<pop(i).fitness&&i~=j
            index=index+1;
        elseif pop(i).fitness==pop(j).fitness&&pop(j).index~=-1&&i~=j
            index=index+1;
        end
    end
    pop(i).index=index;
end
a1=0.6;
for i=1:popsize
    pop(i).fitness=a1*(1-a1)^(pop(i).index-1);
end

function pop=sele(pop)
popsize=size(pop,2);
totalfit=zeros(1,popsize);
for i=1:popsize
    if i==1
       totalfit(i)=pop(i).fitness; 
    else
       totalfit(i)=totalfit(i-1)+pop(i).fitness;
    end
end
totalfit=totalfit/totalfit(popsize);
for i=1:popsize
    p=rand;
    index=1;
    while totalfit(index)<p
        index=index+1;
    end
    new(i)=pop(index);
end
pop=new;


function pop=cr(pop,px,LB,UB)
popsize=size(pop,2);
ndim=length(pop(1).x);
for i=1:popsize
    index=ceil(rand(1,2).*popsize);
    while index(1)==index(2)
        index=ceil(rand(1,2).*popsize);
    end
    if rand<px
       pos=ceil(rand.*ndim);
       v1=pop(index(1)).x(pos);
       v2=pop(index(2)).x(pos);
       pop(index(1)).x(pos)=rand*v2+(1-rand)*v1;
       pop(index(2)).x(pos)=rand*v1+(1-rand)*v2;
       pop(index(1)).x(pos)=boundtest(pop(index(1)).x(pos),LB(pos),UB(pos));
       pop(index(2)).x(pos)=boundtest(pop(index(1)).x(pos),LB(pos),UB(pos));  
    end
end

function pop=mut(pop,pm,LB,UB,iterm,iterm_max)
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

function pop=in(popsize,LB,UB)     %初始化函数
NC=size(LB,1);
for i=1:popsize
    for j=1:NC
        z1=rand;
        while z1==0.25||z1==0.5||z1==0.75
           z1=rand;
        end
        for k=1:200
          z1=4*z1*(1-z1);
        end 
        pop(i).x(j)=LB(j,1)+(UB(j,1)-LB(j,1))*z1;
   end
end




