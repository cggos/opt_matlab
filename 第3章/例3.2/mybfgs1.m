function [xmin,minf]=mybfgs1(fun,hfun,gfun,phi,hphi,gphi,x0,mu,lam,sigma,x_syms,esp)   %DEPÀ„∑®
n=length(x0);
h0=eye(n);
gk=dlag1(hfun,gfun,phi,hphi,gphi,x0,mu,lam,sigma,x_syms);
k=1;
while k<2000
  if norm(gk)<=esp
     xmin=x0;
     minf=fun(xmin);
     break
  end
  d0=-h0\gk;
  m=0;maxm=20;beta=0.5;sigma1=0.2;
  y1=0;
  while m<=maxm
      x1=x0+beta^m*d0';
      x2=lag1(fun,hfun,gfun,x1,mu,lam,sigma); 
      x3=lag1(fun,hfun,gfun,x0,mu,lam,sigma);
      x4=dlag1(hfun,gfun,phi,hphi,gphi,x0,mu,lam,sigma,x_syms);
      x5=x3+sigma1*beta^m*x4'*d0;
      if x2<=x5
          y1=m;  
          break;
      end
      m=m+1;
  end
  y=beta^y1;
  xk=x0'+y*d0;
  gk1=dlag1(hfun,gfun,phi,hphi,gphi,xk',mu,lam,sigma,x_syms);
  sk=xk-x0';
  yk=gk1-gk;
  if (yk'*sk)>0
      h0=h0-(h0*(sk*sk')*h0)/(sk'*h0*sk)+(yk*yk')/(yk'*sk);
  end
  gk=gk1;
  x0=xk';
  k=k+1;
end