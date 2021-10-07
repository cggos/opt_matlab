function y=mydiff1(fun,varargin)
if nargin==1
    n=1;
    syms x
    xsyms=sym(x);
elseif nargin==2
    xsyms=varargin{1};
    n=1;
elseif nargin==3
    xsyms=varargin{1};
    n=varargin{2};
end
num=length(xsyms);
for i=1:num
   y{1}=diff(fun,n,xsyms(i));
end
