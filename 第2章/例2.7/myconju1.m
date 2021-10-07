function [xmin,minf]=myconju1(fun,phi,x0,x_syms,xsyms,esp)   %共轭梯度法求极值
if nargin<6
    esp=1e-6;
end
z=myjacobian1(phi,x_syms);
n=length(x_syms);
g0=eval(subs(z,xsyms,x0));
s0=-g0;
k=0;
while k<5000
   if norm(g0)<=esp
        xmin=x0;
        if ~isempty(fun)
           minf=fun(xmin);
        else
           minf=eval(subs(phi,x_syms,xmin));
        end
        %minf=fun(xmin);
        break
   end
   for i=1:n
      if ~isempty(fun) 
          x3=mysearch1(fun,phi,x0,x_syms,s0,'d');
      else
         x3=mysearch1([],phi,x0,x_syms,s0,'d');
      end
     %x3=mysearch1(fun,phi,x0,x_syms,xsyms,s0,'d');  %一维搜索求步长
     if x3==1
         if ~isempty(fun) 
           x3=mysearch1(fun,phi,x0,x_syms,s0,'equ');
         else
           x3=mysearch1([],phi,x0,x_syms,s0,'equ');
         end
        % x3=mysearch1(fun,phi,x0,x_syms,xsyms,s0,'equ'); 
     end
     xk=x0+x3*s0;
     gk=eval(subs(z,xsyms,xk)); 
     if norm(gk)<=esp
         xmin=xk;
         if ~isempty(fun)
           minf=fun(xmin);
         else
           minf=eval(subs(phi,x_syms,xmin));
         end
         %minf=fun(xmin);
         break
     else
         miuk=(norm(gk))^2/(norm(g0))^2;
         sk=-gk+miuk*s0;
         x0=xk;
         s0=sk;
         g0=gk;
      end
   end
   if exist('xmin','var')==1
        break
   else
       x0=xk;
       g0=eval(subs(z,xsyms,x0));
       s0=-g0; 
   end
   k=k+1;
 end