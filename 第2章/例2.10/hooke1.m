function [xmin,minf]=hooke1(fun,x0,x_syms,d,alpha,esp)  %hooke直接搜索法，x0为初始值，以行输入,d为步长
if nargin==5
    esp=1e-8;
elseif nargin==4
    alpha=1.0;
    esp=1e-8;
elseif nargin==3
    d=1;
    alpha=1.0;
    esp=1e-8;
end
n=size(x0,1);
x(1,:)=x0;
y(1,:)=x(1,:);
h=eye(n);
k=1;
while k<2000
  j=1;
  while 1
     a=y(j,:)+d*h(j,:);
     if isa(fun,'function_handle')
         fa=fun(a);
         fy=fun(y(j,:));
     elseif isa(fun,'sym')
         fa=eval(subs(fun,x_syms,a));
         fy=eval(subs(fun,x_syms,y(j,:)));
     end
     if fa<fy
        y(j+1,:)=a;
     else
       a=y(j,:)-d*h(j,:);
       if isa(fun,'function_handle')
         fa=fun(a);
         fy=fun(y(j,:));
       elseif isa(fun,'sym')
         fa=eval(subs(fun,x_syms,a));
         fy=eval(subs(fun,x_syms,y(j,:)));
       end
       if fa<fy
          y(j+1,:)=a;
       else
          y(j+1,:)=y(j,:);  
       end
     end
     if j<n
         j=j+1;
     else
         break;
     end
  end
  if isa(fun,'function_handle')
     fy=fun(y(n+1,:));
     fx=fun(x(k,:));
  elseif isa(fun,'sym')
     fy=eval(subs(fun,x_syms,y(n+1,:)));
     fx=eval(subs(fun,x_syms,x(k,:)));
  end
  if fy<fx
     x(k+1,:)=y(n+1,:);
     y(1,:)=x(k+1,:)+alpha*(x(k+1,:)-x(k,:));
     k=k+1;
  else
       if d<=esp
          xmin=x(k,:);
          if isa(fun,'function_handle')
             minf=fun(xmin);
          elseif isa(fun,'sym')
             minf=eval(subs(fun,x_syms,xmin));
         end
         break;
       else
          y(1,:)=x(k,:);
          x(k+1,:)=x(k,:);
          d=d/2;
          k=k+1;
        end
  end
end
             
            
        
    
  
