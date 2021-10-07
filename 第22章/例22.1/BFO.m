function [best_x,best_f]=BFO(fun,bacterialnum,step,nc,ns,nre,ned,ped,sr,LB,UB)%细菌觅食算法,求极大
NC=size(LB,1);
for i=1:bacterialnum
    bacterial(i).x=LB'+(UB-LB)'.*rand(1,NC);
    bacterial(i).fitness=fun(bacterial(i).x);
    y(i)=bacterial(i).fitness;
end
[a,b]=sort(y,'descend');
cbest=bacterial(b(1));
bacterial=bacterial(b);
J=y;
n=0;flag=0;
while 1
  for l=1:ned
    for n=1:nre
        for m=1:nc
           for i=1:bacterialnum
               temp_bacterial=bacterial(i);
               delta=unifrnd(-1,1,1,NC);
               bacterial(i).x=bacterial(i).x+step.*delta./norm(delta);
               bacterial(i).x=boundtest(bacterial(i).x,LB,UB);
               bacterial(i).fitness=fun(bacterial(i).x);
               t=0;
               while t<ns
                   t=t+1;
                   if bacterial(i).fitness>=temp_bacterial.fitness
                       temp_bacterial=bacterial(i);
                       bacterial(i).x=bacterial(i).x+step.*delta./norm(delta);
                       bacterial(i).x=boundtest(bacterial(i).x,LB,UB);
                       bacterial(i).fitness=fun(bacterial(i).x);
                   else
                       t=ns;
                       bacterial(i)=temp_bacterial;
                   end
               end
               y(i)=bacterial(i).fitness;
               J(i)=J(i)+y(i);
           end
        end
        [a,b]=sort(J,'descend');     %以能量函数衡量
        %[a,b]=sort(y,'descend');    %以适应度衡量
        bacterial=bacterial(b);
        for i=1:floor(bacterialnum*sr)
            bacterial(bacterialnum-i+1)=bacterial(i);
        end
    end
    for i=1:bacterialnum
        if ped>rand
           bacterial(i).x=LB'+(UB-LB)'.*rand(1,NC);
           bacterial(i).fitness=fun(bacterial(i).x);
           y(i)=bacterial(i).fitness;
        end
    end
  end
  n=n+1;
  [a,b]=max(y);
  if a>cbest.fitness
     cbest=bacterial(b);
  end
  if norm(cbest.x-bacterial(b).x)<1e-8
     best_x=cbest.x;
     best_f=cbest.fitness;
     flag=flag+1;
  else
      flag=0;
  end
  if flag>10||n>20
      break;
  end
end










