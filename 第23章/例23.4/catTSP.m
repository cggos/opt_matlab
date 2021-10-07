function [bestx,bestf]=catTSP(city,catnum,SMP,iter_max)    %猫群算法求TSP
[citynum,d]=city2d(city);
index=randperm(catnum);
for i=1:catnum
    cat(i).route=randperm(citynum);
    if i<=ceil(0.2*catnum) %20%的猫执行跟踪模式
       cat(index(i)).flag=1;
    else
       cat(index(i)).flag=0; 
    end
    point=citynum*rand(1,2);
    while point(1)==point(2)
        point=citynum*rand(1,2);
    end
    v(i,:)=ceil(point);
    y(i)=value(cat(i).route,d);
end
[best_y best_index]=min(y);
z_best=cat(best_index);
y_z_best=best_y;
for i=1:iter_max
    for j=1:catnum
        if cat(j).flag==1   %跟踪行为
           b1=realvelocity(2.0*rand,locationlocation(cat(j).route,z_best.route));
           v1=velocityvelocity(v(j,:),b1);
           cat(j).route=locationvelocity(cat(j).route,v1);
           cat(j).route=isin_TSP(cat(j).route,citynum);
           y(j)=value(cat(j).route,d);
           for k=1:citynum
              n=ceil(citynum*rand);
              while n==k
                 n=ceil(citynum*rand);
              end
              cat1=locationvelocity(cat(j).route,[k n]);
              cat1=isin_TSP(cat1,citynum);
              f1=value(cat(j).route,d);
              if f1<y(j)
                 cat(j).route=cat1;
                 y(j)=f1;
              end
           end
        elseif cat(j).flag==0
           for k=1:SMP
               current_cat(k).route=cat(j).route;
               m1=[];m2=[];
               for t=1:citynum
                  if rand<0.2
                     m1=[m1 k];
                     m2=[m2 current_cat(k).route(t)];
                  end
               end
               if isempty(m1)||length(m1)<1
                   current_cat(k).route=TSPop(current_cat(k).route,'remm');
               elseif length(m1)==2
                  current_cat(k).route(m1(1))=current_cat(k).route(m1(2));
                  current_cat(k).route(m1(2))=current_cat(k).route(m1(1));
               else
                  sr=randperm(length(m1));
                  m2=m2(sr);
                  current_cat(k).route(m1)=m2;
               end
               current_cat(k).route=isin_TSP(current_cat(k).route,citynum);
               y3(k)=value(current_cat(k).route,d);      
           end
           [a,b]=min(y3);
           if a<y(j)
              cat(j).route=current_cat(b).route;
              y(j)=a;
           end
        end
    end
    [a,b]=min(y);
    for j=1:catnum
        if  a<y_z_best
            z_best=cat(b);
            y_z_best=a;
        end
    end
    mr=0.6+0.3*i/iter_max;
    index=randperm(catnum);
    for j=1:catnum
        point=citynum*rand(1,2);
        while point(1)==point(2)
           point=citynum*rand(1,2);
        end
        v(j,:)=ceil(point);
        if j<=ceil(mr*catnum) %2%的猫执行跟踪模式
           cat(index(j)).flag=1;
        else
           cat(index(j)).flag=0; 
        end
    end
end
bestx=z_best.route;
bestf=y_z_best;
TSPplot(city,bestx);


function y=velocityvelocity(v1,v2)    %速度与速度相加
y=[v1;v2];

function y=realvelocity(c,v)    %实数乘以速度
y=floor(c);
n=size(v,1);
if y<1
    y=v(ceil(n*rand),:);
elseif c>size(v,1)
    y=v;
else
    y=v(ceil(n*rand(1,y)),:);
end

function y=locationlocation(x1,x2)   %位置与位置的减法
if length(x1)~=length(x2)
    error('两序列长度不一样');
end
n=length(x1);
if isequal(x1,x2)
    y=ceil(n*rand(1,2));
else
    y=[];
    for i=1:n       
      if x2(i)~=x1(i)
        b=find(x1==x2(i));
        a=[i b];
        y=[y;i b];
        x1=locationvelocity(x1,a);
      end
    end
end

function location=locationvelocity(location,velocity)   %位置与速率的加法
n=size(velocity,1);
for i=1:n
    temp=location(velocity(i,1));
    location(velocity(i,1))=location(velocity(i,2));
    location(velocity(i,2))=temp;
end




    
    
    
    

