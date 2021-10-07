function [xmin,minf]=mygrad1(fun,phi,x0,x_syms,esp)    %梯度法求极值,x_syms为变量
syms lamda
if nargin<6
    esp=1e-6;
end
z=myjacobian1(phi,x_syms);
k=0;
while k<5000  %迭代次数
    f=eval(subs(z,x_syms,x0));
    if norm(f)<=esp
         xmin=x0;
         if ~isempty(fun)
             minf=fun(xmin);
         else
             minf=eval(subs(phi,x_syms,xmin));
         end
         break
    end
    if ~isempty(fun)
        x3=mysearch1(fun,phi,x0,x_syms,-f,'d');%一维搜索求步长 
    else
        x3=mysearch1([],phi,x0,x_syms,-f,'d');
    end
    xk=x0-x3*f;
    x0=xk;
    k=k+1;
end
    
    
        
    