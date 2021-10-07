function [best_x,fval]=SFLA(fun,LB,UB,frogNum,m,max_iterm,L) %混合蛙跳算法，求极小
NC=size(LB,1);
for i=1:frogNum
    frog(i,:)=LB'+(UB-LB)'.*rand(1,NC);
    y(i)=fun(frog(i,:));
end
[a,b]=sort(y,'ascend');  %适应度函数及排序
frog=frog(b,:);
y=y(b);
fval=y(1);best_x=frog(1,:);  %全局最优解，求最小值
for iter=1:max_iterm
  [new,fitness]=grouping(frog,y,m);   %分组
  [new,fitness,best_x,fval]=frog_alter(fun,new,fitness,m,LB,UB,L,best_x,fval);%更新
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

function [new,fitness,best_x,fval]=frog_alter(fun,old,fitness,m,LB,UB,L,best_x,fval)   %更新
NC=length(best_x);
step=(UB-LB)'./100;%每一维的移动步长
for i=1:m     %每个模因组的青蛙
  for j=1:L
     [y_best,b1]=min(fitness(i).y);
     x_best=old(i).x(b1,:); %最优青蛙
     [y_worst,b2]=max(fitness(i).y);
     x_worst=old(i).x(b2,:); %最差青蛙
     for k=1:NC
        if x_best(k)-x_worst(k)>=0
           d_temp(k)=max([rand*(x_best(k)-x_worst(k)),step(k)]);
        else
           d_temp(k)=min([rand*(x_best(k)-x_worst(k)),-step(k)]); 
        end
     end
    temp=x_worst+d_temp;%更新后的位置
    y_temp=fun(temp);
    if y_temp>y_worst
       new(i).x(b2,:)=temp;
       fitness(i).y(b2)=y_temp;
       if y_temp<fval
           best_x=temp;
           fval=y_temp;
       end
    else
       x_best=best_x;
       for k=1:NC
           if x_best(k)-x_worst(k)>=0
              d_temp(k)=max(0.729*rand*2.05*[rand*(x_best(k)-x_worst(k)),step(k)]);
           else
              d_temp(k)=min([0.729*rand*2.05*rand*(x_best(k)-x_worst(k)),-step(k)]); 
           end
       end
       temp=x_worst+d_temp;
       y_temp=fun(temp);
       if y_temp<y_worst
           new(i).x(b2,:)=temp;
           fitness(i).y(b2)=y_temp;
           if y_temp<fval
              best_x=temp;
              fval=y_temp;
           end
       else
           new(i).x(b2,:)=LB'+(UB-LB)'.*rand(1,NC);
           fitness(i).y(b2)=fun(new(i).x(b2,:));
           if fitness(i).y(b2)<fval
               best_x=new(i).x(b2,:);
               fval=fitness(i).y(b2);
           end 
        end
     end
  end
end



