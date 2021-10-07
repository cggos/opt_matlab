function [fun,res]=orthfun(x1,b,m)  %正交多项式最小二乘
syms t
n=length(x1);
h=x1(2)-x1(1);
z=(x1-x1(1))./h;
for i=1:n
  y(i,:)=mylegendre(m,n-1,z(i));
end
for i=1:m+1
    O(:,i)=b(:).*y(:,i);
end
for i=1:m+1
    y1(i)=sum1_1(y(:,i),1);
    y2(i)=sum(O(:,i),1);
end
a=y2./y1;
fun=a(1);
for i=1:m
    [y2,y1]=mylegendre(i,n-1,1);
    num=length(y1);
    a1=1;f=1;
    for j=1:num-1
        a1=a1*(t-j+1);
        f=f+a1*y1(j+1);
    end
    fun=fun+f*a(i+1);
end
fun=simplify(fun);
syms x
fun=subs(fun,t,(x-x1(1))/h);
fun=simplify(fun);
res=0;
for i=1:n
   res=res+(subs(fun,x,x1(i))-b(i))^2;
end
fun=vpa(fun,6);
res=eval(res);


function [y,y1,y2]=mylegendre(m,n,x)     %勒让德多项式
y=ones(1,m+1);
y1=ones(1,m+1);
b1=1;
for i=1:m
    s=1;a=1;b=1;
    for j=1:i
        a=a*(x-j+1);
        b=b*(n-j+1);
        s=s+(-1)^j*factorial(i)/(factorial(j)*factorial(i-j))...
          *factorial(i+j)/(factorial(j)*factorial(i))*a/b;  
    end
    y(1,i+1)=s;
    b1=b1*(n-i+1);
    y1(i+1)=(-1)^i*factorial(m)/(factorial(i)*factorial(m-i))...
          *factorial(m+i)/(factorial(i)*factorial(m))/b;
end
y2=y(end);



function y=sum1_1(x,ndim)
if nargin<2
    ndim=1;
end
[r,c]=size(x);
if r==1
   y=0;
   for i=1:c
      y=y+x(i)^2;
   end
elseif ndim==1  %按列计算
     y=zeros(c,1);
     for j=1:c
        for i=1:r
           y(j)=y(j)+x(i,j)^2;
        end
     end
elseif ndim==2    %按行计算
    y=zeros(r,1);
    for i=1:r
       for j=1:c
          y(i)=y(i)+x(i,j)^2;
       end
    end
end

    

