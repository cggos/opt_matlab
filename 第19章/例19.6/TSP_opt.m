function [new,fitness]=TSP_opt(varargin)  %多步(num）逆转
type=varargin{end};
num=varargin{end-1};
if strcmp(type,'tsp')
    y=varargin{1};
    d=varargin{2}; 
    f=-value(y,d); 
elseif strcmp(type,'bit')
    fun=varargin{1};
    y=varargin{2};
    f=fun(y);
end
new=y;
fitness=f;
k=0;
while 1
    k=k+1;
    if k>num
        break;
    end
    y1=TSPop(y,'em');  %多步逆转
    if strcmp(type,'tsp')
        y2=-value(y1,d);
    else
        y2=fun(y1);
    end
    if y2>f
       new=y1;
       fitness=y2;
       break
    end
end