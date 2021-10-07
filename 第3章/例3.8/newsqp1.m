function [xmin,minf,mu,lam]=newsqp1(fun,hfun,gfun,x0,x_syms,mu,lam)   %SQP求解约束优化
l=size(hfun,1);m=size(gfun,1);
if nargin==5
    mu=zeros(l,1);
    lam=zeros(m,1);
end
n=length(x0);
ro=0.5;eta=0.1;bk=eye(n);sigma=0.8;esp1=1e-6;esp2=1e-5;
hk=eval(subs(hfun,x_syms,x0));
gk=eval(subs(gfun,x_syms,x0));
df=myjacobian1(fun,x_syms);
dfk=eval(subs(df,x_syms,x0))';
dh=myjacobian1(hfun,x_syms);
Ae=eval(subs(dh,x_syms,x0));
dg=myjacobian1(gfun,x_syms);
Ai=eval(subs(dg,x_syms,x0));
k=0;
while k<500
    dx=newquad(bk,dfk,Ae,-hk,Ai,-gk);
    mp1=norm(hk,1)+norm(max(-gk,0),1);
    if norm(dx,1)<esp1&&mp1<esp2
        xmin=x0;
        break;
    end
    deta=0.05;
    tau=max(norm(mu,inf),norm(lam,inf));
    if sigma*(tau+deta)<1
        sigma=sigma;
    else
        sigma=1.0/(tau+2*deta);
    end
    i=0;
    while i<=20
        temp=eta*ro^i*dphi1(fun,hfun,gfun,x_syms,x0,sigma,dx);
        if phi1(fun,hfun,gfun,x_syms,(x0'+ro^i*dx)',sigma)-phi1(fun,hfun,gfun,x_syms,x0,sigma)<temp
            mk=i;
            break
        end
        i=i+1;
        if i==20
            mk=10;
        end
    end
    alpha=ro^mk;
    xk=x0'+alpha*dx;
    hk=eval(subs(hfun,x_syms,xk'));
    gk=eval(subs(gfun,x_syms,xk'));
    dfk=eval(subs(df,x_syms,xk'))';
    Ae=eval(subs(dh,x_syms,xk'));
    Ai=eval(subs(dg,x_syms,xk'));
    Ak=[Ae;Ai];
    lamu=pinv(Ak)'*dfk;
    if l>0&&m>0
        mu=lamu(1:l);
        lam=lamu(l+1:l+m);
    end
    if l==0
        mu=[];lam=lamu;
    end
    if m==0
        mu=lamu;lam=[];
    end
    sk=alpha*dx;
    yk=dlax(fun,hfun,gfun,x_syms,xk',mu,lam)-dlax(fun,hfun,gfun,x_syms,x0,mu,lam);
    if sk'*yk>0.2*sk'*bk*sk
        omega=1;
    else
        omega=0.8*sk'*bk*sk/(sk'*bk*sk-sk'*yk);
    end
    zk=omega*yk+(1-omega)*bk*sk;
    bk=bk+zk*zk'/(sk'*zk)-(bk*sk)*(bk*sk)'/(sk'*bk*sk);
    x0=xk';
    k=k+1;
end
minf=eval(subs(fun,x_syms,xmin));

function y=phi1(fun,hfun,gfun,x_syms,x0,sigma)
f=eval(subs(fun,x_syms,x0));
h=eval(subs(hfun,x_syms,x0));
g=eval(subs(gfun,x_syms,x0));
gn=max(-g,0);
l0=length(h);m0=length(g);
if l0==0
    y=f+1.0/sigma*norm(gn,1);
end
if m0==0
    y=f+1.0/sigma*norm(h,1);
end
if l0>0&&m0>0
   y=f+1.0/sigma*(norm(gn,1)+norm(h,1));
end

function dp=dphi1(fun,hfun,gfun,x_syms,x0,sigma,d)
df1=myjacobian1(fun,x_syms);
df=eval(subs(df1,x_syms,x0))';
h=eval(subs(hfun,x_syms,x0));
g=eval(subs(gfun,x_syms,x0));
gn=max(-g,0);
l0=length(h);m0=length(g);
if l0==0
    dp=df'*d-1.0/sigma*norm(gn,1);
end
if m0==0
    dp=df'*d-1.0/sigma*norm(h,1);
end
if l0>0&&m0>0
    dp=df'-1.0/sigma*(norm(gn,1)+norm(h,1));
end


function y=dlax(fun,hfun,gfun,x_syms,x0,mu,lm)
df1=myjacobian1(fun,x_syms);
df=eval(subs(df1,x_syms,x0))';
dh=myjacobian1(hfun,x_syms);
Ae=eval(subs(dh,x_syms,x0));
dg=myjacobian1(gfun,x_syms);
Ai=eval(subs(dg,x_syms,x0));
[m1,m2]=size(Ai);[l1,l2]=size(Ae);
if l1==0
    y=df-Ai'*lm;
end
if m1==0
    y=df-Ae'*mu;
end
if l1>0&&m1>0
    y=df-Ae'*mu-Ai'*lm;
end





    
    



