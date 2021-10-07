function [xmin,minf]=mylm1(varargin)  %解非线性方程组或非线性最小二乘
if nargin==1    %解方程
   x0=varargin{1};
   esp=1e-6;
   type=1; 
end
if nargin==2    %解方程
   fun=varargin{1};
   x0=varargin{2};
   esp=1e-6;
   type=1; 
end
if nargin==3   %解方程
    fun=varargin{1};
    x0=varargin{2};
    esp=1e-6;
    type=1; 
    x_syms=varargin{3};
end
if nargin==4
    fun=varargin{1};
    x0=varargin{2};
    x_syms=varargin{3};
    esp=varargin{4};
    type=1; 
end
if nargin==6   %最小二乘
    fun=varargin{1};
    x0=varargin{2};
    x_syms=varargin{3};
    esp=1e-6;
    type=2;
    x=varargin{4}; 
    b=varargin{5};
    b_syms=varargin{6};   
end
n=length(x0);
if type==2
   num=size(x,1);
   for i=1:num
       f(i,:)=subs(fun,x_syms,[x(i) b(i)]);    %数据替换
   end
end
muk=norm(Fk(x0));
k=0;
while k<500
        fk=Fk(x0);
        jfk=JFk(x0);
        gk=jfk'*fk;
        dk=-(jfk'*jfk+muk*eye(n))\gk;
        if norm(gk)<esp
            xmin=x0;
            break;
        end
        beta=0.5;sigma=0.4;
        m=0;mk=0;
        while m<20
           fnew=0.5*norm(Fk(x0+beta^m*dk))^2;
           fold=0.5*norm(Fk(x0))^2;
           if fnew<fold+sigma*beta^mk*gk'*dk
              mk=m;
              break;
           end
           m=m+1;
        end
        x0=x0+beta^mk*dk;
        muk=norm(Fk(x0));
        k=k+1;
end
minf=0.5*muk^2;


function F=Fk(x)
n=length(x);
F=zeros(n,1);
for i=1:n
    if i==n
        F(i)=x(1)*x(n)-1;
    else
        F(i)=x(i)*x(i+1)-1;
    end
end


function JF=JFk(x)
n=length(x);JF=zeros(n,n);
for i=1:n-1
    for j=1:n
        if j==i
            JF(i,j)=x(j+1);
        elseif j==i+1
            JF(i,j)=x(j-1);
        end
    end
end
JF(n,1)=x(n);JF(n,n)=x(1);



        