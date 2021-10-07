function cbest=CAT(fun,catnum,SMP,SRD,iter_max,LB,UB)   %猫群算法,求极大
NC=size(LB,1);
cat=myinitialize_cat(catnum,LB,UB);
for i=1:catnum
    cat(i).fitness=fun(cat(i).x);
    y(i)=cat(i).fitness;
end
[a,b]=max(y);  %求最优值
cbest=cat(b);
flag=0;
for iter=1:iter_max
  for i=1:catnum
      if cat(i).flag==1   %跟踪行为
         for j=1:NC
             cat(i).v(j)=cat(i).v(j)+2.0*rand*(cbest.x(j)-cat(i).x(j));
             if cat(i).v(j)>1
                cat(i).v(j)=1;
             elseif cat(i).v(j)<-1
                 cat(i).v(j)=-1;
             end
             cat(i).x(j)=cat(i).x(j)+cat(i).v(j);
             cat(i).x(j)=boundtest(cat(i).x(j),LB(j,1),UB(j,1));
         end
      else
          for k=1:SMP
              current_cat(k)=cat(i);
          end
          for k=1:SMP
              for j=1:NC
                 current_cat(k).x(j)=current_cat(k).x(j)*(1+SRD*(rand*2-1));
              end
              current_cat(k).x=boundtest(current_cat(k).x,LB,UB);
              current(k).fitness=fun(current_cat(k).x);
              y1(k)=current(k).fitness;
          end
          [a,b1]=sort(y1,'descend');
          for k=1:SMP
              p(k)=abs(cat(k).fitness-a(end))/(a(1)-a(end));
          end
          for k=1:SMP
             p1=rand;
             index=1;
             while p(index)<p1
                index=index+1;
                if index>SMP
                    flag=1; 
                    break;
                end
             end
             if flag==0
                cat(i).x=current_cat(index).x;
                cat(i).fitness=current_cat(index).fitness;
             else
                cat(i).x=current_cat(b1(1)).x;
                cat(i).fitness=current_cat(b1(1)).fitness;
             end
          end
      end
      cat(i).fitness=fun(cat(i).x); 
      y(i)=cat(i).fitness;
   end
   [a,b]=max(y);
   if a>cbest.fitness
      cbest=cat(b);
   end
   mr=0.6+0.3*iter/iter_max;
   index=randperm(catnum);
   for i=1:catnum
      if i<=ceil(mr*catnum) %2%的猫执行跟踪模式
         cat(index(i)).flag=1;
      else
         cat(index(i)).flag=0; 
      end
   end
end




function cat=myinitialize_cat(catnum,LB,UB)
NC=size(LB,1);
index=randperm(catnum);
for i=1:catnum
    for j=1:NC
      cat(i).x(j)=LB(j,1)+(UB(j,1)-LB(j,1))*rand;
      cat(i).v(j)=rand;
    end
    if i<=ceil(0.2*catnum) %2%的猫执行跟踪模式
       cat(index(i)).flag=1;
    else
       cat(index(i)).flag=0; 
    end
end







