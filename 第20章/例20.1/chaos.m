function [best_x,fval]=chaos(fun,LB,UB,varargin)   %混沌算法
type=varargin{end};
if type==1   %基本混沌算法
    if nargin==4
        N1=10000;
        N2=500;
    else
        N1=varargin{1};
        N2=varargin{2};
    end
elseif type==2
    if nargin==4
        N1=10000;
        N2=30;
        N3=500;
    else
        N1=varargin{1};
        N2=varargin{2};
        N3=varargin{3};
    end
end
NC=size(LB,1);
fval=inf;
for j=1:NC
    z1=rand;
    while z1==0.25||z1==0.5||z1==0.75
       z1=rand;
    end
    x1(j)=4*z1*(1-z1);
end
for i=1:N1  %粗搜索
   mx=LB'+(UB-LB)'.*x1;%变换到自变量的区间
   y=fun(mx);
   if y<fval
      best_x=mx;
      fval=y;
   end
   x1=4*x1.*(1-x1);  %产生n1维的混沌变量
end
if type==1
   for i=1:N2     %细搜索
      alpha=1-((i-1)/i)^2;
      x=best_x+alpha.*(x1-0.5);
      y=fun(x);
      if y<fval
         best_x=x;
         fval=y;
      end
      x1=4*x1.*(1-x1);
   end
elseif type==2
   gama=unifrnd(0,0.5);
   alpha=0.02;
   for i=1:N2
       for j=1:NC
          LBi(j,1)=best_x(j)-gama*(UB(j,1)-LB(j,1));
          if LBi(j,1)<LB(j,1)
              LBi(j,1)=LB(j,1);
          end
          UBi(j,1)=best_x(j)+gama*(UB(j,1)-LB(j,1));
          if UBi(j,1)<UB(j,1)
              UBi(j,1)=UB(j,1);
          end
       end
       LB=LBi;UB=UBi;
       x=(best_x-LB')./(UB-LB)';
       for k=1:N3
          y=(1-alpha).*x+alpha.*x1;
          z=LB'+(UB-LB)'.*y;
          y1=fun(z);
          if y1<fval
             best_x=z;
             x=(best_x-LB')./(UB-LB)';
             fval=y1;
          end
          x1=4*x1.*(1-x1);
       end
       alpha=alpha/2;
   end
end