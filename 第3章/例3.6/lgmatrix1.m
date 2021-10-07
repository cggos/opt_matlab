function y=lgmatrix1(fun,hfun,x0,mu,x_syms)  %¿≠∏Ò¿ »’æÿ’Û
l=length(mu);
d1=lghess1(fun,hfun,x0,mu,x_syms);
dh=myjacobian1(hfun,x_syms);
d2=eval(subs(dh,x_syms,x0));
y=[d1 -d2';-d2 zeros(l,l)];