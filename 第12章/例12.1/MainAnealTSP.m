function  [best_fval,best_route]=MainAnealTSP(varargin)
if nargin==1
   cityposition=varargin{1};
   flag=0;
else
   cityposition=varargin{1};
   city=varargin{2};       %固定起点城市序号
   flag=1;
end
if size(cityposition,2)>2
   d=cityposition;
else
   d=squareform(pdist(cityposition)); 
end        
n=size(d,1);%the number of cities
route=[];
fval=inf;      %求初始解
rand('state',sum(clock));
for j=1:1000
  if flag==1
     route0=randperm(n);
     a=find(route0==city);
     route0=redu(route0,a,'c');
     route=[city route0];
  else
     route0=randperm(n);
  end
  temp=value(route0,d);
  if temp<fval
     route=route0;
     fval=temp;
  end
end
best_route=route;
best_fval=fval;
t=1;%初始温度
for i=1:30000
     [df,route_after]=exchange(route,d,flag,3);
     if  df<0
         route=route_after;
         fval=value(route,d);
     elseif exp(-df/t)>rand
         route=route_after;
         fval=value(route,d);
     end
     if fval<best_fval
        best_route=route_after;
        best_fval=fval;
     end    
     t=0.999*t;
     if t<0.1^30  %最低温度 
         break
     end
end
if size(cityposition,2)==2     %输入为城市坐标
   TSPplot(cityposition,best_route);
end

