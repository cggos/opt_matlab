function y=greed(varargin)   %贪心算法改善旧路程
if nargin==2
    route=varargin{1};
    d=varargin{2};
else
    route=varargin{1};
    d=varargin{2};
    point=varargin{3};
end
NC=length(route);
y=route;
if nargin==2
   point=[];
   for i=2:NC
      temp=d(y(i-1),:);
      temp(y(1:i-1))=inf;
      [a,y(i)]=min(temp);    %与一个城市最近的城市
   end
else
   temp=d(y(point),:);
   temp(y(point))=inf;
   [a,b]=min(temp);
   y=matrixinsert(y,b,point);
end
