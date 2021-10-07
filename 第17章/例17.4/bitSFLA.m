function [bestx,fval]=bitSFLA(fun,frogNum,m,max_iterm,L,dimension) 
%混合蛙跳算法，求0-1规划,dimension为维数即物件的数目
for i=1:frogNum
    frog(i,:)=-5+10.*rand(1,dimension);
    c(i,:)=realbit(frog(i,:));   %二进制
    y(i)=fun(c(i,:));
end
[a,b]=sort(y,'ascend');  %适应度函数及排序
frog=frog(b,:);y=y(b);c=c(b,:);
fval=y(1);best_x=frog(1,:);  %全局最优解，求最小值
bestx=c(1,:);
for iter=1:max_iterm
  [new,fitness]=grouping(frog,y,m);   %分组
  [new,fitness,best_x,fval,bestx]=frog_alter2(fun,new,fitness,m,L,best_x,fval,bestx);%更新
  old=[];
  y=[];
  for i=1:m
      old=[old;new(i).x];   %重组
      y=[y fitness(i).y];
  end
  [a,b]=sort(y,'ascend');  %适应度函数及排序
  frog=old(b,:);
  y=y(b);
end

function [new,fitness]=grouping(old,y,m)   %分组
N=size(old,1);
for i=1:m    %分组
    j=1;
    while 1
        if m*(j-1)+i>N
            break
        else
           new(i).x(j,:)=old(m*(j-1)+i,:);   
           fitness(i).y(j)=y(m*(j-1)+i);
           j=j+1;
        end
    end
end

function [old,fitness,best_x,fval,bestx]=frog_alter2(fun,old,fitness,m,L,best_x,fval,bestx)   %更新
NC=length(best_x);
for i=1:m     %每个模因组的青蛙
  for j=1:L
     [y_best,b1]=min(fitness(i).y);
     x_best=old(i).x(b1,:); %最优青蛙
     [y_worst,b2]=max(fitness(i).y);
     x_worst=old(i).x(b2,:); %最差青蛙
     for k=1:NC
        d_temp(k)=rand*(x_best(k)-x_worst(k));
     end
    temp=x_worst+d_temp;%更新后的位置
    temp1=realbit(temp);
    y_temp=fun(temp1);
    if y_temp<y_worst
       old(i).x(b2,:)=temp;
       fitness(i).y(b2)=y_temp;
       if y_temp<fval
           best_x=temp;
           fval=y_temp;
           bestx=temp1;
       end
    else
       x_best=best_x;
       for k=1:NC
          d_temp(k)=0.729*rand*2.05*rand*(x_best(k)-x_worst(k));
       end
       temp=x_worst+d_temp;
       temp1=realbit(temp);
       y_temp=fun(temp1);
       if y_temp<y_worst
           old(i).x(b2,:)=temp;
           fitness(i).y(b2)=y_temp;
           if y_temp<fval
              best_x=temp;
              fval=y_temp;
              bestx=temp1;
           end
       else
           old(i).x(b2,:)=-5+10.*rand(1,NC);
           fitness(i).y(b2)=fun(old(i).x(b2,:));
           if fitness(i).y(b2)<fval
               best_x=old(i).x(b2,:);
               fval=fitness(i).y(b2);
               bestx=realbit(best_x);
           end 
        end
     end
  end
end

function b=realbit(x)   %实数转化成二进制
n=length(x);
b=zeros(1,n);
for i=1:n
    if 1/(1+exp(-x(i)))<=0.5
        b(i)=0;
    else
        b(i)=1;
    end
end



