function [xmin,minf]=mysucfail1(fun,x0,h,esp)   % �ɹ���ʧ�ܷ���ֵ,hΪ��������
if nargin<4
    esp=1e-6;
end
f0=fun(x0);
n=0;
while abs(h)>esp
   n=n+1;
   if n>5000
        error('�޽�');
   end
   x1=x0+h;
   if isa(fun,'function_handle')
     f1=fun(x1);
   elseif isa(fun,'sym')
     f1=eval(subs(fun,x1));
   end   
 %  f1=fun(x1);
   if f1<f0
     x0=x1;
     f0=f1;
     h=2*h;
   else
     h=-h/4;
   end
end
xmin=x0;
if isa(fun,'function_handle')
   minf=fun(xmin);
elseif isa(fun,'sym')
   minf=eval(subs(fun,xmin));
end   
%minf=fun(xmin);