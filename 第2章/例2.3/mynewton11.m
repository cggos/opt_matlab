function [xmin,minf]=mynewton11(fun,phi,varargin)  %牛顿法求极值,fun为函数计算公式，phi为符号
if nargin==3
   x0=varargin{1};
   esp=1e-6;
end
if nargin==4
   x0=varargin{1};
   esp=varargin{2};
end
ds1=mydiff1(phi);
ds2=mydiff1(ds1);
%ds1=char(ds1);
%ds2=char(ds2);
x1=x0-eval(subs(ds1,x0))/eval(subs(ds2,x0));
n=1;
while abs(x1-x0)>esp
    n=n+1;
    if n>5000
       error('无解');
    end
    x0=x1;
    x1=x0-eval(subs(ds1,x1))/eval(subs(ds2,x1)); 
end
xmin=x0;
minf=fun(xmin);