function [bestx,fval]=gradnewton(fun,fun1,iter_max,x0,x_syms,esp)   %最速下降法－修正牛顿法
if nargin<6
    esp=1e-6;
end
Z=myjacobian1(fun1,x_syms);
n=length(x0);
for i=1:iter_max
    k=0;
    while k<100  %迭代次数
        f=eval(subs(fun,x_syms,x0));
        if norm(f)<=esp
            xmin=x0;
            minf=f;
            break
        end
        for j=1:n
            if x0(j)==0
               Dx(j)=0.00001;
            else
               Dx(j)=0.00001*x0(j);
            end
        end
        total=0;
        for j=1:n
            x0(j)=x0(j)+Dx(j);
            DF(j)=((eval(subs(fun,x_syms,x0))-f)/Dx(j));
            total=total+DF(j)*DF(j);
            x0(j)=x0(j)-Dx(j);
        end
        RL=f/total;
        for j=1:n
           x0(j)=x0(j)-RL*DF(j);
        end
        k=k+1;
    end
    if exist('xmin','var')==0
       F1=eval(subs(fun,x_syms,x0));
    else
        F1=minf;
    end
    new=x0;
    k=0;
    while k<500
        f=eval(subs(fun1,x_syms,x0));
        z=eval(subs(Z,x_syms,x0));
        if norm(f)<=esp
           break
        end
        for j=1:5
            x0=x0'-inv(z)*f;
            x0=x0';
            f=eval(subs(fun1,x_syms,x0));   
        end
        k=k+1;
    end
    minf=eval(subs(fun,x_syms,x0));
    if minf>=F1
       bestx=new;
       fval=F1;
       break
    end
end

