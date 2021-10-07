function cbest=myga(fun,numvar,popsize,iterm_max,pm,px,LB,UB)   %遗传算法
if isempty(pm)
    pm=0.1;
end
if isempty(px)
    px=0.9;
end
if isempty(iterm_max)
    iterm_max=5000;
end
if isempty(popsize)
    popsize=30;
end
pop=in(numvar,popsize,LB,UB);
for i=1:popsize
    pop(i).fitness=fun(pop(i).x);
end
cbest=pop(1);
cworst=pop(1);
[cbest,cworst]=f(pop,1,cbest,cworst);   %最优与最差
for iterm=1:iterm_max
    pop=ca(pop);    %适应度值
    pop=s(pop);     %选择
    pop=c(pop,px,LB,UB); %交叉
    pop=m(pop,pm,LB,UB,iterm,iterm_max); %变异
    for i=1:popsize
        pop(i).fitness=fun(pop(i).x);
    end
    [cbest,cworst,best]=f(pop,iterm,cbest,cworst);
    pop(popsize)=best;
end


function pop=in(numvar,popsize,LB,UB)     %初始化函数，L为每个基因串的长度
for i=1:popsize 
    pop(i).x=LB'+rand(1,numvar).*(UB-LB)';
end
 
 function [cbest,cworst,best]=f(pop,iterm,cbest,cworst)
popsize=size(pop,2);
best=pop(1);
worst=pop(1);
for i=2:popsize
    if pop(i).fitness<best.fitness
        best=pop(i);
    elseif pop(i).fitness>worst.fitness
        worst=pop(i);
    end
end
if iterm==1
    cbest=best;
    cbest.index=1;
else
    if best.fitness<cbest.fitness
        cbest=best;
        cbest.index=iterm;
    end
end

function pop=ca(pop)
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

function pop=s(pop)
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

function pop=c(pop,px,LB,UB)
popsize=size(pop,2);
ndim=length(pop(1).x);
for i=1:popsize
    index=ceil(rand(1,2).*popsize);
    if rand<px
       pos=ceil(rand.*ndim);
       v1=pop(index(1)).x(pos);
       v2=pop(index(2)).x(pos);
       pop(index(1)).x(pos)=rand*v2+(1-rand)*v1;
       pop(index(2)).x(pos)=rand*v1+(1-rand)*v2;
       pop(index(1)).x(pos)=b(pop(index(1)).x(pos),LB(pos),UB(pos));
       pop(index(2)).x(pos)=b(pop(index(1)).x(pos),LB(pos),UB(pos));  
    end
end

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

function y=b(x,LB,UB)    %边界检测
if x>=UB
    y=LB+(x-UB);
elseif x<=LB
    y=UB-(LB-x);
else
    y=x;
end
    

