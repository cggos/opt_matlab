function y=BB(fun,hfun,x0,mu,tau,x_syms)    %SQP方法中的求KKT矩阵
h=lghess1(fun,hfun,x0,mu,x_syms);
dhfun=myjacobian1(hfun,x_syms);
dh=eval(subs(dhfun,x_syms,x0));
y=h+1.0/(2*tau)*(dh'*dh);