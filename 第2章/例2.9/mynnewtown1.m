function  [xmin,minf]=mynnewtown1(fun,phi,x0,x_syms,esp)    %ÄâÅ£¶Ù·¨
if nargin<4
    esp=1e-6;
end
n=length(x0);
hk=eye(n);
z=myjacobian1(phi,x_syms);
k=1;
while k<2000
    gk=eval(subs(z,x_syms,x0))';
    if norm(gk)<=esp
         xmin=x0;
         if ~isempty(fun)
             minf=fun(xmin);
         else
             minf=eval(subs(phi,x_syms,xmin));
         end
         %minf=fun(xmin);
         break
    end
    dk=-hk*gk;
    if ~isempty(fun)
        x3=mysearch1(fun,phi,x0,x_syms,dk','a');
    else
        x3=mysearch1([],phi,x0,x_syms,dk','a');
    end
    %x3=mysearch1(fun,phi,x0,x_syms,xsyms,dk','a');  
    xk=x0'+x3*dk;
    sk=xk-x0';
    yk=eval(subs(z,x_syms,xk'))'-gk;
    hk=hk+(sk-hk*yk)*(sk-hk*yk)'/((sk-hk*yk)'*yk);
    x0=xk';
    k=k+1;
end

    