function [best_x,fval]=qGA3(fun)   %改进量子遗传算法,求极大值
%num为每个基因段的长度
prompt={'量子数';'最大迭代数';'变量下界';'变量上界';'变异概率';'自变量离散精度'};
name='输入算法各参数';
defaultanswer={'30','200','-inf','inf','0.1','1e-8'};
answer=inputdlg(prompt,name,1,defaultanswer);
popsize=str2num(answer{1});
iter_max=str2num(answer{2});
LB=str2num(answer{3});
UB=str2num(answer{4});
pm=str2num(answer{5});
eps=str2num(answer{6});
NC=size(LB,1);
%num=NC*max(ceil(log2((UB-LB)/eps + 1)));      %根据自变量离散精度，确定二进制编码位串的长度
num=max(ceil(log2((UB-LB)/eps + 1)));  
for i=1:popsize
    pop(i).p=zeros(1,NC*num);
    pop(i).fai=pi/4*ones(1,NC*num);
    pop(i).q=zeros(2,NC*num);
    for j=1:NC*num
        pop(i).q(1,j)=cos(pop(i).fai(j));
        pop(i).q(2,j)=sin(pop(i).fai(j)); 
    end
    pop(i).fitness=0;
end
for i=1:popsize
    for j=1:NC*num
        p=rand;
        if p>0.5
            pop(i).p(1,j)=1;
        else
            pop(i).p(1,j)=0;
        end
    end
end
pop=Dec1(pop,LB,UB,NC*num,1);%解码,1为变量域的序号
for i=1:popsize
   pop(i).fitness=fun(pop(i).var);
   f(i)=pop(i).fitness;
end
[a,b]=max(f);   %极大值
cbest=pop(b);
table=[0.0001*pi  -1  1  round(rand*2-1) 0;
       0.0001*pi  1  -1  round(rand*2-1) 0;
       0.0001*pi  1  -1  0               1;
       0.0001*pi  1  -1  round(rand*2-1) 0;
       0.0001*pi -1  1   round(rand*2-1) 0;
      0.0001*pi  1  -1   0               round(rand*2-1);
      0.0001*pi -1   1   0               round(rand*2-1);
      0.0001*pi -1   1   0               round(rand*2-1)];
flag=0;   %控制变量
for iter=2:iter_max 
  pop=crossover_ABC_TSP(pop,NC,num);
  pop=Dec1(pop,LB,UB,NC*num,1);%解码,1为变量域的序号
  for i=1:popsize
     pop(i).fitness=fun(pop(i).var);
     f(i)=pop(i).fitness;
  end
  [a,b]=max(f);   %极大值
  cbest1=pop(b);
  favg=mean(f);
  if cbest1.fitness>cbest.fitness
      cbest=cbest1;
  end
  for i=1:popsize
     if flag==3 
         delta=0.001*pi+((cbest.fitness-pop(i).fitness)/cbest.fitness+exp(1-iter_max/iter)/2)*0.049*pi;
         flag=0;
     else
        delta=0.001*pi+((cbest.fitness-pop(i).fitness)/cbest.fitness)*0.049*pi;
     end
    % if pop(i).fitness>=favg
     %    delta=0.001*pi*(favg-pop(i).fitness)/(cbest.fitness-favg);
     %else
     %    delta=0.05*pi;
     %end
     delta_sita=0;
     for j=1:NC*num
       s=0;
       if pop(i).p(j)==0&&cbest.p(j)==0
         if pop(i).fitness<cbest.fitness
            delta_sita=table(1,1);
         else
            delta_sita=table(2,1);
         end
       elseif pop(i).p(j)==0&&cbest.p(j)==1
         if pop(i).fitness<cbest.fitness
            delta_sita=table(3,1);
         else
           delta_sita=table(4,1);
           if pop(i).q(1,j)*pop(i).q(1,j)>0
              s=table(4,2);
           elseif  pop(i).q(1,j)*pop(i).q(1,j)<0
              s=table(4,3);
           elseif  pop(i).q(1,j)==0
              s=table(4,4);
           elseif  pop(i).q(2,j)==0
              s=table(4,5);
           end
         end
       elseif pop(i).p(j)==1&&cbest.p(j)==0
         if pop(i).fitness<cbest.fitness
            delta_sita=table(5,1);
            if pop(i).q(1,j)*pop(i).q(1,j)>0
              s=table(5,2);
            elseif  pop(i).q(1,j)*pop(i).q(1,j)<0
              s=table(5,3);
            elseif  pop(i).q(1,j)==0
              s=table(5,4);
            elseif  pop(i).q(2,j)==0
              s=table(5,5);
            end
          else
            delta_sita=table(6,1);
            if pop(i).q(1,j)*pop(i).q(1,j)>0
               s=table(6,2);
            elseif  pop(i).q(1,j)*pop(i).q(1,j)<0
              s=table(6,3);
            elseif  pop(i).q(1,j)==0
              s=table(6,4);
            elseif  pop(i).q(2,j)==0
              s=table(6,5);
            end
          end
        elseif pop(i).p(j)==1&&cbest.p(j)==1
          if pop(i).fitness<cbest.fitness
            delta_sita=table(7,1);
            if pop(i).q(1,j)*pop(i).q(1,j)>0
              s=table(7,2);
            elseif  pop(i).q(1,j)*pop(i).q(1,j)<0
              s=table(7,3);
            elseif  pop(i).q(1,j)==0
              s=table(7,4);
            elseif  pop(i).q(2,j)==0
              s=table(7,5);
            end
           else
             delta_sita=table(8,1);
             if pop(i).q(1,j)*pop(i).q(1,j)>0
               s=table(8,2);
             elseif  pop(i).q(1,j)*pop(i).q(1,j)<0
               s=table(8,3);
             elseif  pop(i).q(1,j)==0
               s=table(8,4);
             elseif  pop(i).q(2,j)==0
               s=table(8,5);
             end
           end
       end
       newpop=pop(i);
       %newpop.fai(j)=pop(i).fai(j)+s*delta_sita;
       newpop.fai(j)=pop(i).fai(j)+s*delta;
       if rand<pm
          newpop.fai(j)=newpop.fai(j)+pi/2;
       end
       newpop.q(1,j)=cos(newpop.fai(j));
       newpop.q(2,j)=sin(newpop.fai(j));
     end
       pop(i).q=newpop.q; 
       for j=1:NC*num
          p=rand;
          if p>pop(i).q(j)*pop(i).q(j)
            pop(i).p(1,j)=1;
          else
            pop(i).p(1,j)=0;
          end
        end
  end
  pop=Dec1(pop,LB,UB,NC*num,1);%解码
  for k=1:popsize
      pop(k).fitness=fun(pop(k).var);
      f(k)=pop(k).fitness;
  end
  [a,b]=max(f);
  cbest1=pop(b);
  if cbest1.fitness>cbest.fitness
     cbest=cbest1;
     flag=0;
  else
      flag=flag+1;
  end
end
fval=cbest.fitness;
best_x=cbest.var;

function pop=crossover_ABC_TSP(pop,NC,num)  %量子相干交叉算子
popsize=size(pop,2);
for i=1:popsize
    start=1;
    fin=num;
    k=randperm(popsize);
    for j=1:NC
      pop(i).p(start:fin)=pop(k(j)).p(start:fin);
      pop(i).fai(start:fin)=pop(k(j)).fai(start:fin);
      pop(i).q(1,start:fin)=pop(k(j)).q(1,start:fin);
      pop(i).q(2,start:fin)=pop(k(j)).q(2,start:fin);
      start=start+num;
      fin=fin+num;
    end
end