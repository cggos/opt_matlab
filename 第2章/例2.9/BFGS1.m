function [xmin,minf]=BFGS1(fun,phi,x0,x_syms,esp)   %BFGSÀ„∑®
n=length(x0);
h0=eye(n);
z=myjacobian1(phi,x_syms);
gk=eval(subs(z,x_syms,x0))';
k=1;
while k<2000
  if norm(gk)<=esp
     xmin=x0;
     if ~isempty(fun)
        minf=fun(xmin);
     else
         minf=eval(subs(phi,x_syms,xmin));
     end
     break
  end
  d0=-h0\gk;
  if ~isempty(fun)
     lamda=mysearch1(fun,phi,x0,x_syms,d0','a');
  else
     lamda=mysearch1([],phi,x0,x_syms,d0','a');
  end
  xk=x0'+lamda*d0;
  gk1=eval(subs(z,x_syms,xk'))';
  sk=xk-x0';
  yk=gk1-gk;
  if (yk'*sk)>0
      h0=h0-(h0*(sk*sk')*h0)/(sk'*h0*sk)+(yk*yk')/(yk'*sk);
  end
  gk=gk1;
  x0=xk';
  k=k+1;
end