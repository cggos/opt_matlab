function [best_x,fval]=bitABC(varargin)%离散蜂群算法
type=varargin{end};
if strcmp(type,'tsp')   %求TSP
    beenum=varargin{1};
    iter_max=varargin{2};
    limit=varargin{3};
    city=varargin{4};
    [NC,d]=city2d(city);
elseif strcmp(type,'bit')    %求0－1规划或背包问题
    fun=varargin{1};
    beenum=varargin{2};
    iter_max=varargin{3};
    limit=varargin{4};  
    NC=varargin{5};
end
x=init_ABC(beenum,NC,type);
for i=1:beenum
    bee(i).x=x(i,:);
    if strcmp(type,'tsp') 
       y(i)=-value(x(i,:),d); 
    elseif strcmp(type,'bit')
       y(i)=fun(x(i,:));
    end
    bee(i).oBas=0;
    bee(i).fitness=y(i);
    f(i)=y(i);
end
[a,b]=max(f);
gbest=bee(b);
[new,fitness]=grouping(x,y,2);
for j=1:size(new(1).x,1)
   employbee(j).x=new(1).x(j,:);    %引领蜂(雇佣蜂)
   employbee(j).oBas=0;
   employbee(j).fitness=fitness(1).y(j);
end
for j=1:size(new(2).x,1)
   onlookbee(j).x=new(2).x(j,:);    %跟随蜂(观察蜂）
   onlookbee(j).oBas=0;
   onlookbee(j).fitness=fitness(2).y(j);
end
pbest=onlookbee;
for iter=2:iter_max
    for i=1:length(employbee)
      if employbee(i).oBas>limit    %引领蜂
          if strcmp(type,'tsp')
             y5=randperm(NC);
             y5=greed(y5,d);
          elseif strcmp(type,'bit')
             y5=rand(1,NC)<0.5; 
          end
          y5=TSPop(y5,'em');   %逆转算子
          y5=TSPop(y5,'imm');  %免疫算子
          if strcmp(type,'tsp')    %多步逆转
             [new1,fitness]=TSP_opt(y5,d,10,'tsp');
          elseif strcmp(type,'bit')
              [new1,fitness]=TSP_opt(fun,y5,10,'bit');
          end
          employbee(i).x=new1;
          employbee(i).oBas=0;
          employbee(i).fitness=fitness;
      else
          neighbour=ceil(rand*beenum);   %邻域
          while 1
              if isequal(bee(neighbour).x,employbee(i).x)||neighbour>beenum
                 neighbour=ceil(rand*beenum);
              else
                  break;
              end
          end
          y6(1,:)=TSPop(employbee(i).x,bee(neighbour).x,'cpmx');   %部分匹配交叉算子
          y7=TSPop(y6(1,:),'em');   %逆转算子
          y7=TSPop(y7,'imm');  %免疫算子
          if strcmp(type,'tsp')  %多步逆转 
             [new1,fitness]=TSP_opt(y7,d,1,'tsp');   
          elseif strcmp(type,'bit')
              [new1,fitness]=TSP_opt(fun,y7,1,'bit');
          end
          df=employbee(i).fitness-fitness;
          if df<0
              employbee(i).x=new1;
              employbee(i).oBas=0;
              employbee(i).fitness=fitness;
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
          if strcmp(type,'tsp')
             y5=randperm(NC);
          elseif strcmp(type,'bit')
             y5=rand(1,NC)<0.5; 
          end
          y5=TSPop(y5,'em');   %逆转算子
          y5=TSPop(y5,'imm');  %免疫算子
          if strcmp(type,'tsp')    %多步逆转
             [new1,fitness]=TSP_opt(y5,d,5,'tsp');
          elseif strcmp(type,'bit')
              [new1,fitness]=TSP_opt(fun,y5,5,'bit');
          end
          onlookbee(i).x=new1;
          onlookbee(i).oBas=0;
          onlookbee(i).fitness=fitness;
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
          y6(1,:)=TSPop(onlookbee(i).x,employbee(neighbour).x,'cpmx');   %部分匹配交叉算子
          y7=TSPop(y6(1,:),'em');   %逆转算子
          y7=TSPop(y7,'imm');  %免疫算子
          if strcmp(type,'tsp')    %多步逆转
             [new1,fitness]=TSP_opt(y7,d,5,'tsp');
          elseif strcmp(type,'bit')
              [new1,fitness]=TSP_opt(fun,y7,5,type);
          end
          df=onlookbee(i).fitness-fitness;
          if df<0
              onlookbee(i).x=new1;
              onlookbee(i).oBas=0;
              onlookbee(i).fitness=fitness;
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
    if a>gbest.fitness
        gbest.fitness=a;
        gbest.x=x2(b,:);
    end   
end
best_x=gbest.x;
fval=gbest.fitness;
if strcmp(type,'tsp')
  TSPplot(city,best_x);
end


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

