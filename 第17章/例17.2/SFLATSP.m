function [best_x,fval]=SFLATSP(city,frogNum,m,max_iterm,L,pc) %混合蛙跳算法
[NC,d]=city2d(city);
for i=1:frogNum
    frog(i,:)=randperm(NC);
    y(i)=value(frog(i,:),d);
end
[a,b]=sort(y,'ascend');  %适应度函数及排序
frog=frog(b,:);
y=y(b);
fval=y(1);best_x=frog(1,:);  %全局最优解
for iter=2:max_iterm
  [new,fitness]=grouping(frog,y,m);   %分组
  new=frog_alterTSP(new,fitness,NC,m,L,best_x,d);%更新
  old=[];
  y=[];
  for i=1:m
      old=[old;new(i).x];   %重组
  end
  old=mutation_TSP(old,d,pc);
  for i=1:frogNum
     y(i)=value(old(i,:),d);
  end
  [a,b]=sort(y,'ascend');  %适应度函数及排序
  frog=old(b,:);
  y=y(b);
  fval=y(1);best_x=frog(1,:);
end
TSPplot(city,best_x)

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

function old=frog_alterTSP(old,fitness,NC,m,L,best_x,d)
for i=1:m
  for j=1:L
     [y_best,b1]=min(fitness(i).y);
     x_best=old(i).x(b1,:); %最优青蛙
     [y_worst,b2]=max(fitness(i).y);
     x_worst=old(i).x(b2,:); %最差青蛙
     for k=1:NC
        if x_best(k)-x_worst(k)>=0
           temp(k)=x_worst(k)+min([ceil(rand*(x_best(k)-x_worst(k))),3]);%更新后的位置
        else
           temp(k)=x_worst(k)+max([ceil(rand*(x_best(k)-x_worst(k))),-3]); 
        end
     end
     temp=isin_TSP(temp);%检验是否有重复或缺少
     y_temp=value(temp,d);
     if y_temp<y_worst
         old(i).x(b2,:)=temp;
    else
       x_best=best_x;
       for k=1:NC
         if x_best(k)-x_worst(k)>=0
           temp(k)=x_worst(k)+min([ceil(rand*(x_best(k)-x_worst(k))),3]);
         else
           temp(k)=x_worst(k)+max([ceil(rand*(x_best(k)-x_worst(k))),-3]); 
         end
       end
       temp=isin_TSP(temp);
       y_temp=value(temp,d);
       if y_temp<y_worst
          old(i).x(b2,:)=temp;
       else
          old(i).x(b2,:)=randperm(NC);
       end
     end
  end
end

function y=isin_TSP(y)   %将重复的城市划去，补上缺的城市
NC=length(y);
lost=[];
for i=1:NC
    if isin(i,y)==0           %缺的城市
      lost=[lost i];
    end
end
m=find(y>NC);    %超出范围
if ~isempty(m)
   m=[m find(y<=0)];
else
   m=find(y<=0);
end
if ~isempty(m)
    for i=1:length(m)
        y(m(i))=lost(i);   
    end
    lost=redu(lost,1:length(m),'c');
end
if ~isempty(lost)
    for i=1:NC
       a=find(y==i);
       if length(a)>1
           for j=1:length(a)-1
               y(a(j))=lost(j);
           end
           lost=redu(lost,1:length(a)-1,'c');
       end
       if isempty(lost)
           break
       end
    end
end

function y=isin(x,A)   %判断是否在A中
k=0;
for i=1:length(A)
   if abs(x-A(i))<=0.0001
      k=k+1;
      y(k)=i;
      break
   else
       y=0;
   end
end

function route=mutation_TSP(route,d,pc)
[d_min,d_index]=sort(d,2);%按列排列
d_index=d_index(:,2);%与各城市最近的城市
[num,NC]=size(route);
for i=1:num
    temp1=route(i,:);
    y_temp1=value(temp1,d);
    for j=1:NC
      p=rand;
      if p<pc
          a=route(i,j);   %要变异的城市
          temp=d_index(a,1);%与a最近的城市
          b=find(route(i,:)==temp);
          if j==NC
             aa=route(i,1);
             route(i,1)=route(i,b);
             route(i,b)=aa;
          else
             aa=route(i,j+1);
             route(i,j+1)=route(i,b);
             route(i,b)=aa;
          end      
      end
    end
    temp2=route(i,:);
    y_temp2=value(temp2,d);
    if y_temp2>y_temp1
       route(i,:)=temp1;
    end
end






 

 

