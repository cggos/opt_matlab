function [x_max,y_max]=ABC(fun,varargin)%蜂群算法
beenum=varargin{1};
iter_max=varargin{2};
LB=varargin{3};
UB=varargin{4};
limit=varargin{5};
NC=size(LB,1);
for i=1:beenum
    bee(i).x=LB'+(UB-LB)'.*rand(1,NC);
    bee(i).oBas=0;
    bee(i).fitness=fun(bee(i).x);
end
[y_max,x_max]=my_max(bee);
for iter=1:iter_max
    for i=1:beenum
      if bee(i).oBas>limit
           bee(i).x=LB'+(UB-LB)'.*rand(1,NC);
           bee(i).oBas=0;
           bee(i).fitness=fun(bee(i).x);
      end     
    end
    y=f1_ABC(bee);%选择概率　
    followbeenum=ceil(beenum.*y);
    newbee=bee(1);
    for i=1:beenum
       for j=1:followbeenum(i)
           neighbour=ceil(rand*beenum);
           while neighbour==i
               neighbour=ceil(rand*beenum);
           end
           newbee.x=bee(i).x+abs(bee(neighbour).x-bee(i).x).*unifrnd(-1,1,1,NC);
           newbee.x=boundtest(newbee.x,LB,UB);
           newbee.fitness=fun(newbee.x);
           if newbee.fitness<bee(i).fitness
               bee(i).oBas=bee(i).oBas+1;
           else
               bee(i).oBas=0;
               bee(i)=newbee;
           end
       end
    end
   [y_max1,x_max1]=my_max(bee);
   if y_max1>y_max
      y_max=y_max1;
      x_max=x_max1;
   end
end

function [y_max,x_max]=my_max(bee)
num=size(bee,2);
for i=1:num
   y1(i)=bee(i).fitness;
end
[y_max,b]=max(y1);
x_max=bee(b).x;

function y=f1_ABC(bee)
num=size(bee,2);
y=zeros(1,num);
a=0;
for i=1:num
    a=a+bee(i).fitness;
end
for i=1:num
    y(i)=bee(i).fitness/a;
end