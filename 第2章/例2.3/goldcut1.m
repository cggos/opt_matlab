function [x0,minf]=goldcut1(fun,varargin)   %黄金分割法求最小点
if nargin==4
    a=varargin{1};
    b=varargin{2};
    esp=varargin{3};
end
if nargin==3
    a=varargin{1};
    b=varargin{2};
    esp=1e-6;
end
if nargin==2
    esp=1e-6;
    [a,b]=interval1(fun,varargin{1},0.1,'f');
end
n=0;
while (abs((a-b))>esp)
    if n>5000
       error('无解');
    end
    n=n+1;
    x1=a+(b-a)*0.618;
    if isa(fun,'function_handle')
       y1=fun(x1);
    elseif isa(fun,'sym')
       y1=eval(subs(fun,x1));
    end   
    %y1=fun(x1);
    x2=b-(b-a)*0.618;
    if isa(fun,'function_handle')
       y2=fun(x2);
    elseif isa(fun,'sym')
       y2=eval(subs(fun,x2));
    end   
    %y2=fun(x2);
    if y1<y2
        a=x2;
    elseif y1>y2
        b=x1;
    elseif y1==y2
        a=x2;
        b=x1;
    end  
end
x0=(a+b)/2;
if isa(fun,'function_handle')
   minf=fun(x0);
elseif isa(fun,'sym')
   minf=eval(subs(fun,x0));
end   
%minf=fun(x0);
