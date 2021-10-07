function [x0,minf]=myDF1(phi,varargin)  %对分法求极值
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
    [a,b]=interval1(phi,varargin{1},0.1,'f');
end
L0=b-a;
n=0;
y=2^(-n)*L0+(1-2^(-n))*esp;   %区间长
while y>esp
   if n>5000
      error('无解');
  end
  n=n+1;
  x1=a+(L0-esp)/2;
  x2=b-(L0-esp)/2;
  f1=phi(x1);
  f2=phi(x2);
  L0=(L0+esp)/2;
  if f1<f2
      b=x2;
  elseif f1>f2
      a=x1;
  elseif  f1==f2
      a=x1;
      b=x2;
  end
  y=2^(-n)*L0+(1-2^(-n))*esp;
end
x0=(a+b)/2;
minf=phi(x0);


