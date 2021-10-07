function [bestx,bestf]=mypsoTSP(city)    %粒子群算法求TSP
prompt={'粒子数';'最大迭代数';};
name='输入算法各参数';
defaultanswer={'30','1000'};
answer=inputdlg(prompt,name,1,defaultanswer);
popsize=str2num(answer{1});
max_iterm=str2num(answer{2});
[citynum,d]=city2d(city);
for i=1:popsize
    pop(i,:)=randperm(citynum);
    point=citynum*rand(1,2);
    while point(1)==point(2)
        point=citynum*rand(1,2);
    end
    v(i,:)=ceil(point);
    y(i)=value(pop(i,:),d);
end
[best_y best_index]=min(y);
z_best=pop(best_index,:);
g_best=pop;
y_g_best=y;
y_z_best=best_y;
w_max=0.9;w_min=0.4;
w_min1=0.2;
for i=1:max_iterm
    w=w_max-i*(w_max-w_min)/max_iterm;
    for j=1:popsize
        a1=realvelocity(w,v(j,:));
        b1=realvelocity(0.2*rand,locationlocation(pop(j,:),g_best(j,:)));
        b2=realvelocity(0.2*rand,locationlocation(pop(j,:),z_best));
        v1=velocityvelocity(a1,b1);
        v2=velocityvelocity(v1,b2);
        pop(j,:)=locationvelocity(pop(j,:),v2);
        y(j)=value(pop(j,:),d);
        for k=1:citynum
           n=ceil(citynum*rand);
           while n==k
               n=ceil(citynum*rand);
           end
           pop1=locationvelocity(pop(j,:),[k n]);
           f1=value(pop(j,:),d);
           if f1<y(j)
               pop(j,:)=pop1;
               y(j)=value(pop(j,:),d);
           end
        end
    end
    for j=1:popsize
        if y(j)<y_g_best(j)
            g_best(j,:)=pop(j,:);
            y_g_best(j)=y(j);
        end
        if y(j)<y_z_best
            z_best=pop(j,:);
            y_z_best=y(j);
        end
    end
end
bestx=z_best;bestf=y_z_best;
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




    
    
    
    

