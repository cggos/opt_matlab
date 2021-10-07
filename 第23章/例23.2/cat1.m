function x_max=cat1(fun,catnum,SMP,SRD,iter_max,LB,UB)   %猫群算法,一次求方程的根
NC=size(LB,1);
cat=myinitialize_cat(catnum,LB,UB);
y=cat_fit(fun,cat);
[y_max,x_max]=cat_max(cat,y);%求最优值
for iter=1:iter_max
   for i=1:catnum
      if cat(i).flag==1
         for j=1:NC
             x_best=x_max(j);
             cat(i).v(j)=cat(i).v(j)+2.5*rand*(x_best-cat(i).v(j));
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
          for k=1:SMP-1
             for j=1:NC
                current_cat(k).x(j)=current_cat(k).x(j)+SRD*(rand*2-1);
                current_cat(k).x(j)=boundtest(current_cat(k).x(j),LB(j,1),UB(j,1));               
             end
          end
          y1=cat_fit(fun,current_cat);
          [y_max1,x_max1]=cat_max(current_cat,y1);
          for j=1:NC
              cat(i).x(j)=x_max1(j);
          end 
      end
   end
   y=cat_fit(fun,cat);
   [y_max2,x_max2]=cat_max(cat,y);
   for j=1:NC
     if y_max2(j)>y_max(j)
        y_max(j)=y_max2(j);
        x_max(j)=x_max2(j);
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
    if i<=ceil(0.02*catnum) %2%的猫执行跟踪模式
       cat(index(i)).flag=1;
    else
       cat(index(i)).flag=0; 
    end
end


function y=cat_fit(fun,cat)%一次求解7个根
NC=length(cat(1).x);
r=size(cat,2);
for i=1:r
    for j=1:NC
       y(i,j)=fun(cat(i).x(j));
       y(i,j)=1/(1+y(i,j)^2);%适应度函数
    end
end

function [y_max,x_max]=cat_max(cat,y)
NC=length(cat(1).x);
[y_max,b]=max(y);   %因为是平行求解,所以互相不干扰
for i=1:NC
  x_max(i)=cat(b(i)).x(i);
end











