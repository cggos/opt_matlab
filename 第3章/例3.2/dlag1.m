function dy=dlag1(hfun,gfun,phi,hphi,gphi,x0,mu,lam,sigma,x_syms) 
dfun=myjacobian1(phi,x_syms);
dhfun=myjacobian1(hphi,x_syms);
dgfun=myjacobian1(gphi,x_syms);
he=hfun(x0);
gi=gfun(x0);
dy=eval(subs(dfun,x_syms,x0))';  
dhe=eval(subs(dhfun,x_syms,x0))';
dgi=eval(subs(dgfun,x_syms,x0))';
l=length(he);m=length(gi); 
if ~isempty(dhe)
    for i=1:l
       dy=dy+(sigma*he(i)-mu(i))*dhe(:,i);
    end
else
     dy=dy;
end
if ~isempty(dgi)
   for i=1:m
     dy=dy+(sigma*gi(i)-lam(i))*dgi(:,i);
   end
else
    dy=dy;
end
