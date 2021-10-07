function [x0,minf]=myparabola1(fun,x,esp)   %抛物线法求极值
if nargin<3
    esp=1e-6;
end
if length(x)==1
    [x1,x3]=interval1(fun,x,0.01,'f');
    x2=(x1+x3)/2;
elseif length(x)==2
    x1=x(1);x3=x(2);
    x2=(x1+x3)/2;
else
    x1=min(x);
    x2=max(x);
    x3=median(x);
end
if isa(fun,'function_handle')
  f1=fun(x1);
  f2=fun(x2);
  f3=fun(x3);
elseif isa(fun,'sym')
  f1=eval(subs(fun,x1));
  f2=eval(subs(fun,x2));
  f3=eval(subs(fun,x3));
end
n=0;
while 1
    n=n+1;
    if n>5000
        error('无解');
    end
   c1=(x2^2-x3^2)*f1+(x3^2-x1^2)*f2+(x1^2-x2^2)*f3;
   c2=(x2-x3)*f1+(x3-x1)*f2+(x1-x2)*f3;
   if c2==0
       x0=x2;
       if isa(fun,'function_handle')
          minf=fun(x0);
       elseif isa(fun,'sym')
          minf=eval(subs(fun,x0));
       end   
       break
   end
   xp=0.5*c1/c2;
   if isa(fun,'function_handle')
      fp=fun(xp);
   elseif isa(fun,'sym')
      fp=eval(subs(fun,xp));
   end   
   if abs(x2-xp)<=esp
     if abs(f2-fp)<=esp
        if fp<=f2
            x0=xp;
            if isa(fun,'function_handle')
               minf=fun(x0);
            elseif isa(fun,'sym')
               minf=eval(subs(fun,x0));
            end   
            break;
        else
            x0=x2;
            if isa(fun,'function_handle')
               minf=fun(x0);
            elseif isa(fun,'sym')
               minf=eval(subs(fun,x0));
            end   
            break;
        end
     end
   else
     if fp<=f2
       if xp<=x2
           x3=x2;
           x2=xp;
           f3=f2;
           f2=fp;
       else
           x1=x2;
           x2=xp;
           f1=f2;
           f2=fp;
       end
   else
       if xp<=x2
          x1=xp;
          f1=fp;
       else
          x3=xp;
          f3=fp;
       end
     end
   end
end


   
   
   