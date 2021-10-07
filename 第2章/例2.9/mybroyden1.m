function [xmin,minf]=mybroyden1(fun,phi,x0,x_syms,esp)   %broydenËã·¨
if nargin<4
    esp=1e-6;
end
n=length(x0);
hk=eye(n);ph=0.5;
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
         break
    end
    dk=-hk*gk;
    if ~isempty(fun)
       x3=mysearch1(fun,phi,x0,x_syms,dk','a');
    else
       x3=mysearch1([],phi,x0,x_syms,dk','a'); 
    end
    xk=x0'+x3*dk;
    sk=xk-x0';
    yk=eval(subs(z,x_syms,xk'))'-gk;
    hy=hk*yk;
    sy=sk'*yk;
    yhy=yk'*hk*yk;
    if sy<0.2*yhy
        theta=0.8*yhy/(yhy-sy);
        sk=theta*sk+(1-theta)*hy;
        sy=0.2*yhy;
    end
    vk=sqrt(yhy)*(sk/sy-hy/yhy);
    hk=hk-(hy*hy')/yhy+(sk*sk')/sy+ph*(vk*vk');
    x0=xk';
    k=k+1;
end