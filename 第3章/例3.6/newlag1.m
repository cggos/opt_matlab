function [xmin,minf]=newlag1(fun,hfun,x0,x_syms,esp)  %牛顿－拉格朗日法
if nargin==4
    esp=1e-6;
end
n=length(x0);
l=length(hfun);    %乘子向量的长度
mu=0.2*ones(1,l);
beta=0.5;sigma=0.2;
k=0;
while k<500
    x1=dlg11(fun,hfun,x0,mu,x_syms);   %拉格朗日函数的梯度
    if norm(x1)<esp
        xmin=x0;
        minf=eval(subs(fun,x_syms,x0));
        break;
    end
    x2=lgmatrix1(fun,hfun,x0,mu,x_syms);   %拉格朗日矩阵
    dz=(-x2\x1)';
    dx=dz(1:n);du=dz(n+1:n+l);
    m=0;mk=0;
    while m<20
        t1=beta^m;
        if (norm(dlg11(fun,hfun,x0+t1*dx,mu+t1*du,x_syms))^2<=(1-sigma*t1)*norm(x1)^2)
            mk=m;
            break;
        end
        m=m+1;
    end
    xk=x0+beta^mk*dx;
    mu=mu+beta^mk*du;
    x0=xk;
    k=k+1;
end



    