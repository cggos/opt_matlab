function [route,f]=qGATSP(city)%量子遗传算法求解TSP
[NC,d]=city2d(city);
prompt={'量子数';'最大迭代数';'变异概率'};
name='输入算法各参数';
defaultanswer={'20','2000','0.05'};
answer=inputdlg(prompt,name,1,defaultanswer);
popsize=str2num(answer{1});
iter_max=str2num(answer{2});
pm=str2num(answer{3});
num=length(dec2bin(NC));
for i=1:popsize
    pop(i).p=zeros(1,NC*num);
    pop(i).fai=pi/4*ones(1,NC*num);
    pop(i).q=zeros(2,NC*num);
    pop(i).fitness=0;
end
for i=1:popsize
    for j=1:NC*num
        pop(i).q(1,j)=cos(pop(i).fai(j));
        pop(i).q(2,j)=sin(pop(i).fai(j)); 
    end
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
pop=decode_QGA_TSP(pop,NC,num);%解码
for i=1:popsize
    y(i)=value(pop(i).var,d);
    pop(i).fitness=y(i);
end
[a,b]=sort(y,'ascend');
pop=pop(b);
cbest=pop(b(1));
table=[0 0 0 0 0;...
       0 0 0 0 0;...
       0 0 0 0 0;...
       0.05*pi -1 1 round(rand*2-1) 0;
       0.01*pi -1 1 round(rand*2-1) 0;
       0.025*pi 1 -1 0 round(rand*2-1);
       0.005*pi 1 -1 0 round(rand*2-1);
       0.025*pi 1 -1 0 round(rand*2-1)];
for iter=2:iter_max
    pop=selection_ABC_TSP(pop);    %选择算子
    pop=crossover_ABC_TSP(pop,NC,num);    %交叉算子
    %pop=crossover_qGA(pop,num);
    pop=decode_QGA_TSP(pop,NC,num);    %解码
    for i=1:popsize
       y(i)=value(pop(i).var,d);
       pop(i).fitness=y(i);
    end
    [a,b]=sort(y,'ascend');
    pop=pop(b);
    cbest1=pop(1);
    if cbest1.fitness<cbest.fitness
       cbest=cbest1;
    end
    for i=1:popsize
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
            newpop.fai(j)=pop(i).fai(j)+s*delta_sita;
           if rand<pm      %变异算子
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
    pop=decode_QGA_TSP(pop,NC,num);%解码
    for i=1:popsize
        y(i)=value(pop(i).var,d); 
        pop(i).fitness=y(i);
    end
    [a,b]=sort(y,'ascend');
    pop=pop(b);
    cbest1=pop(1);
    if cbest1.fitness<cbest.fitness
       cbest=cbest1;
    end
end
route=cbest.var;
f=value(route,d);
TSPplot(city,route);

function pop=decode_QGA_TSP(pop,NC,num)%解码，将二位码变成十进制数
popsize=size(pop,2);
for i=1:popsize
    start=1;
    fin=num;
    for j=1:NC
        tvars(1:num)=pop(i).p(start:fin);
        start=start+num;
        fin=fin+num;
        temp1=num2str(tvars);
        temp1=bin2dec(temp1);
        if temp1>NC
            temp1=ceil(rand*NC);
        end
        if temp1<1
            temp1=ceil(rand*NC);
        end
        pop(i).var(j)=temp1;
    end
   pop(i).var=isin_TSP(pop(i).var);
end


function pop=selection_ABC_TSP(pop)
popsize=size(pop,2);
cfitness=zeros(1,popsize);
for i=1:popsize
    if i==1
        cfitness(i)=pop(i).fitness;
    else
        cfitness(i)=cfitness(i-1)+pop(i).fitness;
    end
end
cfitness=cfitness./cfitness(popsize);
for i=1:popsize
    p=rand;
    index=1;
    while cfitness(index)<p
        index=index+1;
    end
    newpop(i)=pop(index);
end
pop=newpop;

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


  






            

