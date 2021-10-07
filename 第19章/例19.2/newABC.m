function [best_x,fval]=newABC(fun,varargin)%蜂群算法新程序
beenum=varargin{1};
iter_max=varargin{2};
LB=varargin{3};
UB=varargin{4};
limit=varargin{5};
NC=size(LB,1);
for i=1:beenum
    for j=1:NC
        x1(i,j)=LB(j,1)+(UB(j,1)-LB(j,1))*rand;
        x2(i,j)=LB(j,1)+UB(j,1)-x1(i,j);   %反向点
    end
    y1(i,1)=fun(x1(i,:));
    y2(i,1)=fun(x2(i,:));
end
y=[y1;y2];
x=[x1;x2];
[a,b]=sort(y,'descend');
y=y(b(1:beenum));
x=x(b(1:beenum),:);   %选最优
for i=1:beenum
   bee(i).x=x(i,:);
   bee(i).oBas=0;
   bee(i).fitness=y(i);
end
%gbest=bee(b(1));
gbest=bee(1);
[new,fitness]=grouping(x,y,2);
for j=1:size(new(1).x,1)
   employbee(j).x=new(1).x(j,:);    %引领蜂
   employbee(j).oBas=0;
   employbee(j).fitness=fitness(1).y(j);
end
for j=1:size(new(2).x,1)
   onlookbee(j).x=new(2).x(j,:);    %跟随蜂
   onlookbee(j).oBas=0;
   onlookbee(j).fitness=fitness(2).y(j);
end
pbest=onlookbee;
for iter=2:iter_max
    for i=1:length(employbee)
      if employbee(i).oBas>limit    %引领蜂
          employbee(i).x=LB'+(UB-LB)'.*rand(1,NC);
          employbee(i).oBas=0;
          employbee(i).fitness=fun(employbee(i).x);
      else
          neighbour=ceil(rand*beenum);   %邻域
          while 1
              if isequal(bee(neighbour).x,employbee(i).x)||neighbour>beenum
                 neighbour=ceil(rand*beenum);
              else
                  break;
              end
          end
          newemploybee.x=gbest.x+abs(bee(neighbour).x-employbee(i).x).*unifrnd(-1,1,1,NC);
          newemploybee.x=boundtest(newemploybee.x,LB,UB,2);
          newemploybee.fitness=fun(newemploybee.x);
          df=employbee(i).fitness-newemploybee.fitness;
          if df<0
              employbee(i).x=newemploybee.x;
              employbee(i).fitness=newemploybee.fitness;
              employbee(i).oBas=0;
          else
              employbee(i).oBas=employbee(i).oBas+1;
          end
      end
    end
    p=f1_ABC(employbee);%选择概率
    [a,b]=sort(p,'descend');
    p=p(b);
    employbee=employbee(b);
    for i=1:length(onlookbee)
       if onlookbee(i).oBas>limit    %观察蜂
          onlookbee(i).x=LB'+(UB-LB)'.*rand(1,NC);
          onlookbee(i).oBas=0;
          onlookbee(i).fitness=fun(onlookbee(i).x);
      else
          p1=rand;
          index=1;
          flag=1;
          while p(index)<p1
              index=index+1;
              if index>length(employbee);
                  flag=0;
                  break;
              end   
          end
          if flag==1
             neighbour=index;   %邻域
          else
             neighbour=ceil(length(employbee)/3*rand);
          end
          newonlookbee.x=pbest(i).x+abs(employbee(neighbour).x-onlookbee(i).x).*unifrnd(-1,1,1,NC);
          newonlookbee.x=boundtest(newonlookbee.x,LB,UB,2);
          newonlookbee.fitness=fun(newonlookbee.x);
          df=onlookbee(i).fitness-newonlookbee.fitness;
          if df<0
              onlookbee(i).x=newonlookbee.x;
              onlookbee(i).fitness=newonlookbee.fitness;
              onlookbee(i).oBas=0;
          else
              onlookbee(i).oBas=onlookbee(i).oBas+1;
          end
       end
       if pbest(i).fitness<onlookbee(i).fitness
           pbest(i)=onlookbee(i);
       end
    end
    for i=1:length(employbee)
        bee(i)=employbee(i);
        x2(i,:)=employbee(i).x;
        y2(i,1)=employbee(i).fitness;
    end
    for i=1:length(onlookbee)
        bee(length(employbee)+i)=onlookbee(i);
        x2(length(employbee)+i,:)=onlookbee(i).x;
        y2(length(employbee)+i,1)=onlookbee(i).fitness;
    end
    [a,b]=max(y2);
    if a>gbest.fitness;
        gbest.fitness=a;
        gbest.x=x2(b,:);
    end
end
best_x=gbest.x;
fval=gbest.fitness;

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