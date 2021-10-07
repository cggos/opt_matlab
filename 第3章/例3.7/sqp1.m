function [xmin,minf]=sqp1(fun,hfun,x0,x_syms,esp)   %SQP方法求解等式约束
if nargin==4
    esp=1e-6;
end
l=size(hfun,1);    %乘子向量的长度
mu=0.2*ones(1,l);
beta=0.6;sigma=0.2;tau=1.55;
k=0;
while k<500
    p1=dlg11(fun,hfun,x0,mu,x_syms);
    if norm(p1)^2<esp
        xmin=x0;
        minf=eval(subs(fun,x_syms,x0));
        break;
    end
    H=BB(fun,hfun,x0,mu,tau,x_syms);
    c=myjacobian1(fun,x_syms);
    c=eval(subs(c,x_syms,x0))';
    Ae=myjacobian1(hfun,x_syms);
    Ae=eval(subs(Ae,x_syms,x0));
    be=-eval(subs(hfun,x_syms,x0));
    [dx,lam]=subf(H,c,Ae,be);
    du=lam-mu'-1.0/(2*tau)*Ae*dx;
    m=0;mk=0;
    while m<20
        p2=dlg11(fun,hfun,x0+beta^m*dx',mu+beta^m*du',x_syms);
        if norm(p2)^2<=(1-sigma*beta^m)*norm(p1)^2
            mk=m;
            break
        end
        m=m+1;
    end
    xk=x0+beta^mk*dx';
    mu=mu+beta^mk*du';
    if norm(xk-x0)<esp
        xmin=x0;
        minf=eval(subs(fun,x_syms,x0));
        break;
    end
    x0=xk;
    k=k+1;
end

