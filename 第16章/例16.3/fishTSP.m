function [best_route,fval]=fishTSP(city,fishnum,max_iterm,visual,try_number,delta)   %人工鱼群算法求TSP
[citynum,d]=city2d(city);
parameter.visual=visual;
parameter.try_number=try_number;
parameter.delta=delta;
for i=1:fishnum
   afish(i,:)=randperm(citynum);
   y(i)=value(afish(i,:),d);
end
[best_val best_num]=min(y);
best_route=afish(best_num,:);      %最优鱼的路径
fval=best_val;
for j=1:max_iterm
    a=exp(-30*(j/max_iterm)^5);
    visual1=parameter.visual*a+0.001;
    parameter.delta=0.95*parameter.delta;
    for i=1:fishnum
       if ~isequal(afish(i,:),best_route)
           parameter.visual=fishdstc(afish(i,:),best_route);
       else
          parameter.visual=visual1; 
       end
       parameter.step=0.90*parameter.visual;
       afish(i,:)=fishevaluate2(afish(i,:),afish,parameter,d);
       y(i)=value(afish(i,:),d);
       if y(i)<fval
           best_route=afish(i,:);
           fval=y(i);
        end
    end
end
TSPplot(city,best_route);


function afish=fishevaluate2(afish0,afish1,parameter,d)
af_follow=fishfollow2(afish0,afish1,parameter,d);
af_swarm=fishswarm2(afish0,afish1,parameter,d);
af_prey=fishprey2(afish0,parameter,d);
afish2=af_follow;
if value(af_swarm,d)<value(afish2,d)
    afish2=af_swarm;
end
if value(af_prey,d)<value(afish2,d)
    afish2=af_prey;
end
if value(afish2,d)<value(afish0,d)
    afish=afish2;
else
    afish=fishmove2(afish0,d);
end

function afish=fishfollow2(afish0,afish1,parameter,d)
fishnum=size(afish1,1);
n=0;
f_min=inf;
min_i=1;
for i=1:fishnum
    if (fishdstc2(afish1(i,:),afish0)<parameter.visual)
        n=n+1;    
       if value(afish1(i,:),d)<f_min
           f_min=value(afish1(i,:),d);
           min_i=i;
       end
    end
end
if (f_min*n<parameter.delta*value(afish0,d))&&(~isequal(afish1(min_i,:),afish0))
    afish=afish1(min_i,:);
else
    afish=fishprey2(afish0,parameter,d);
end

function afish=fishswarm2(afish0,afish1,parameter,d)
fishnum=size(afish1,1);
NC=length(afish0);
n=0;
route_center=1:NC;
for i=1:fishnum
    if (fishdstc2(afish1(i,:),afish0)<parameter.visual)
        n=n+1;
        for j=1:NC
           route_center(j)=route_center(j)+afish1(i,j);
        end
    end
end
if n~=0
   for j=1:NC
       route_center(j)=round(route_center(j)/n);
   end
   center=afish1(1,:);
   for j=2:fishnum
       out1=0;
       out2=0;
       for i=1:NC
           out1=out1+sign(center(i)-route_center(i));
           out2=out2+sign(afish1(j,i)-route_center(i));
       end
       if  out2<out1
           center=afish1(j,:);
       end
   end
   if (value(center,d)*n>value(afish0,d)*parameter.delta)&&(all(center~=afish0))
       afish=center;
       return
   end
end
afish=fishmove2(afish0,d);

function afish=fishprey2(afish0,parameter,d)
NC=length(afish0);
for i=1:parameter.try_number
   loc=sort(unidrnd(NC,1,2));
   while loc(1)==loc(2)
       loc=sort(unidrnd(NC,1,2));
   end
   middle_route=fliplr(afish0(loc(1):loc(2)));
   route_after=[afish0(1:loc(1)-1) middle_route afish0(loc(2)+1:end)];
   if ((value(afish0,d)>value(route_after,d)))
       afish=route_after;
       return
   end
end
afish=fishmove2(afish0,d);

function afish=fishmove2(afish0,d)
NC=length(afish0);
city(1)=1;
i=1;
while i<NC+1
    dmin=inf;
    for k=1:10
        j=ceil(NC*rand);
        while isin(j,city)
          j=ceil(NC*rand);
        end
        if d(i,j)<dmin
           dmin=d(i,j);
           nextcity=j;
        end
    end
    if ~isin(nextcity,city)
       afish(i)=afish0(nextcity);
       city(i)=nextcity;
    else
       i=i-1;
    end
    i=i+1;
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

function y=fishdstc2(afish1,afish2)
NC=length(afish1);
y=0;
for i=1:NC+1
    if i==NC+1
        y=y+sign(afish1(1)-afish2(1));
    else
        y=y+sign(afish1(i)-afish2(i));
    end
end
