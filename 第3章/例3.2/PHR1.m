function [xmin,minf]=PHR1(fun,hfun,gfun,phi,hphi,gphi,x0,x_syms,esp)    %PHRËã·¨£¨³Ë×Ó·¨£©
if nargin==4
    esp=1e-06;
end
sigma=2.0;theta=0.8;eta=2.0;
he=hfun(x0);
gi=gfun(x0);
l=length(he);m=length(gi);
mu=0.1*ones(l,1);lam=0.1*ones(1,m);
betak=10;betaold=10;k=0;
while betak>esp&&k<300
    xmin=mybfgs1(fun,hfun,gfun,phi,hphi,gphi,x0,mu,lam,sigma,x_syms,esp);
    he=hfun(xmin);
    gi=gfun(xmin)';
    betak=sqrt(norm(he,2)^2+norm(min(gi,lam/sigma),2)^2);
    if betak>esp
        if ~isempty(he)
           mu=mu-sigma*he;
        end
        if ~isempty(gi)
           lam=max(0,lam-sigma*gi);
        else
           lam=max(0,lam);
        end
        if (k>=2&&betak>theta*betaold)
            sigma=eta*sigma;
        end
    end
    k=k+1;
    betaold=betak;
    x0=xmin;
end
minf=fun(xmin);



