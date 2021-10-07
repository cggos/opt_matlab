function [pop,cbest,f_avg,cworst]=rotationq(fun,pop,cbest,table,pm,LB,UB)     %量子旋转门
NC=size(LB,1);      %维数
num=length(cbest.x)/NC;    %每个二进制变量长度
popsize=size(pop,2);
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
pop=Dec1(pop,LB,UB,num*NC,1);%解码
for k=1:popsize
    pop(k).fitness=fun(pop(k).var);
    f(k)=pop(k).fitness;
end
[a,b]=sort(f,'ascend');   %从小到大
cbest1=pop(b(1));
cworst=b(end);    %最差个体的序号
if cbest1.fitness<cbest.fitness
   cbest=cbest1;
end
f_avg=mean(f);

        