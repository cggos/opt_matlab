function [xmin,minf]=DFP1(fun,phi,x0,x_syms,esp)   %DEPÀ„∑®
if nargin==5
    esp=1e-6;
end
n=length(x0);
h0=eye(n);
z=myjacobian1(phi,x_syms);
gk=eval(subs(z,x_syms,x0));
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
  d0=-h0*gk';
  if ~isempty(fun)
     lamda=mysearch1(fun,phi,x0,x_syms,d0','a');
  else
      lamda=mysearch1([],phi,x0,x_syms,d0','a');
  end
  xk=x0'+lamda*d0;
  gk1=eval(subs(z,x_syms,xk'));
  dalta1=xk-x0';
  gama1=(gk1-gk)';
  h1=h0+dalta1*dalta1'/(dalta1'*gama1)-h0*gama1*(h0*gama1)'/(gama1'*h0*gama1);
  gk=gk1;
  h0=h1;
  x0=xk';
  k=k+1;
end
  