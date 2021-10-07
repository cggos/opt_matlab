function [minx,minf]=powell1(fun,phi,x0,x_syms,esp)   %powell算法求极值
if nargin==4
    esp=1e-6;
end
num=length(x0);
d=eye(num);   %搜索方向
for k=1:2000
    x1=x0;xk=x0;
    for i=1:num
       if ~isempty(fun)
          y=mysearch1(fun,phi,x0,x_syms,d(i,:),'equ');
       else
           y=mysearch1([],phi,x0,x_syms,d(i,:),'equ');
       end
       x1=x1+y*d(i,:);
       xk=[xk;x1];
    end
    err=norm(xk(end,:)-x0);  
    if err<esp
        minx=xk(end,:);    %最终值
        if ~isempty(fun)
           minf=fun(minx);
        else
            minf=eval(subs(phi,x_syms,minx));
        end
        break;
    else
        if ~isempty(fun)
            y=opfun9(xk);
        else
            n=size(xk,1);
            for i=1:n
               y(i)=eval(subs(phi,x_syms,xk(i,:)));
            end
        end
        for j=1:size(xk,1)-1
           a(j)=y(j)-y(j+1);
        end
        delta=max(a);
        f1=y(1);
        f2=y(end);
        if ~isempty(fun)
            f3=fun(2*xk(end,:)-x0);
        else
            f3=eval(subs(phi,x_syms,2*xk(end,:)-x0));
        end
        if 2*delta<(f1-2*f2+f3)
           x0=xk(end,:);
        else
          d(num,:)=xk(end,:)-x0;
          if ~isempty(fun)
              y1=mysearch1(fun,phi,xk(end,:),x_syms,d(num,:),'equ');
          else
              y1=mysearch1([],phi,xk(end,:),x_syms,d(num,:),'equ');
          end
          x0=xk(end,:)+y1*d(num,1);
        end
    end
end

    
        
    
