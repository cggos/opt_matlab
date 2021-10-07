function [a,b]=myJT1(phi,x0,lamda)  
%进退法求区间,xk为初始点，x1为搜索初始点，lamda搜索步长
if isa(phi,'sym')
    phai0=eval(subs(phi,x0));
else
    phai0=phi(x0);
end
%phai0=phi(x0);
t2=x0+lamda;
if isa(phi,'sym')
    phai2=eval(subs(phi,t2));
else
    phai2=phi(t2);
end
%phai2=phi(t2);
if phai2>phai0
    lamda=-lamda;
    t1=x0+lamda;
    if isa(phi,'sym')
       phai1=eval(subs(phi,t1));
    else
       phai1=phi(t1);
    end
   %phai1=phi(t1);
else
   t1=x0+lamda;
   if isa(phi,'sym')
       phai1=eval(subs(phi,t1));
    else
       phai1=phi(t1);
    end
 %  phai1=phi(t1);
end       
while phai1<=phai0
   lamda=2*lamda;
   t2=x0;
   x0=t1;
   phai0=phai1;
   t1=x0+lamda;
   if isa(phi,'sym')
       phai1=eval(subs(phi,t1));
    else
       phai1=phi(t1);
    end
   %phai1=phi(t1);
end
a=min(t1,t2);
b=max(t1,t2);    

