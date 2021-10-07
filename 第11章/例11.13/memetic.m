function cbest=memetic(fun,popsize,searchnum,pc,pm,iter_max,LB,UB,eps1)%memetic算法
%myval=[80 3 0.6 0.05 100 8];
%popsize=myval(1);searchnum=myval(2);pc=myval(3);pm=myval(4);iter_max=myval(5);num=myval(6);
%LB=[-1;-1];UB=[2;2];eps1=1e-3;
if nargin==8
    eps1=1e-4;
end
NC=size(LB,1);
CodeLen=NC*max(ceil(log2((UB-LB)/eps1 + 1)));  %根据自变量离散精度，确定二进制编码位串的长度 
for i=1:popsize
    pop(i).x=(rand(1,NC*CodeLen)<0.5);     %二进制编码
    pop(i).fitness=0;
end
pop=Dec1(pop,LB,UB,CodeLen);
for i=1:popsize
   pop(i).fitness=fun(pop(i).var);
   f(i)=pop(i).fitness;
end
[a,b]=sort(f);
cbest=pop(b(1));
for i=1:iter_max
    pop=crossover_MTC(pop,pc);
    pop=local_search(fun,pop,searchnum,LB,UB,CodeLen);
    pop=mutation_MTC(pop,pm);
    pop=local_search(fun,pop,searchnum,LB,UB,CodeLen);
    pop=Dec1(pop,LB,UB,CodeLen);
    for j=1:popsize
        pop(j).fitness=fun(pop(j).var);
        f(j)=pop(j).fitness;
    end
    [a,b]=sort(f);
    pop=selection_MTC(pop);
    pop(b(end))=pop(b(1));
    cbest1=pop(b(1));
    if cbest1.fitness<cbest.fitness     %求极小值
         cbest=cbest1;
    end
end

function pop=crossover_MTC(pop,pc)  %交叉算子
popsize=size(pop,2);
NC=length(pop(1).x);
k=randperm(popsize);
for i=1:popsize/2
    p=rand;
    if p<pc
      point=ceil(rand*NC);
      while point==NC
          point=ceil(rand*NC);
      end
      for j=point:NC
          temp=pop(k(2*i-1)).x(j);
          pop(k(2*i-1)).x(j)=pop(k(2*i)).x(j);
          pop(k(2*i)).x(j)=temp;
      end
    end
end

function pop=local_search(fun,pop,searchnum,LB,UB,CodeLen)%memetic算法的局部搜索算子
NC=size(LB,1);
popsize=size(pop,2);
for i=1:popsize
   pop_search=pop(i);
   for t=1:searchnum      
       k=randperm(NC*CodeLen);
       temp=pop_search.x(k(1));
       temp1=pop_search.x(k(2));
       j=2;
       while temp==temp1
           j=j+1;
           temp1=pop_search.x(k(j));   
       end
       pop_search.x(k(1))=pop_search.x(k(j));
       pop_search.x(k(j))=temp;
       pop_search=Dec1(pop_search,LB,UB,CodeLen);
       pop_search.fitness=fun(pop_search.var);
       if pop_search.fitness<pop(i).fitness
            pop(i)=pop_search;
       end
   end
end

function pop=mutation_MTC(pop,pm) %变异算子
NC=length(pop(1).x);
popsize=size(pop,2);
for i=1:popsize
     for j=1:NC
        test=rand;
        if test<pm
            pop(i).x(j)=abs(pop(i).x(j)-1);
        end
    end
end

function pop=Dec1(pop,LB,UB,CodeLen)           %二进制转换为十进制编码
NC=size(LB,1);
popsize=size(pop,2);
sublen=CodeLen/NC;
pow_two=2.^(0:sublen)';
maxintval=((2^sublen))-1;
range=UB'-LB';
for i=1:popsize
   start=1;
   fin=sublen;
   for j=1:NC
      tvars(1:sublen)=pop(i).x(start:fin);
      start=start+sublen;
      fin=fin+sublen;
      temp1=0;
      for k=1:sublen
         temp1=temp1+pow_two(k)*tvars(sublen-k+1);
      end
      pop(i).var(j)=(range(j)*(temp1/maxintval))+LB(j);
   end
end









