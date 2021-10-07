function y=dlg11(fun,hfun,x0,mu,x_syms)  %拉格朗日函数的梯度
dfun=myjacobian1(fun,x_syms);
df=eval(subs(dfun,x_syms,x0));
dhfun=myjacobian1(hfun,x_syms);
h=eval(subs(hfun,x_syms,x0));
num=length(mu);
s=0;
for i=1:num
    s=s+mu(i)*eval(subs(dhfun(1,:),x_syms,x0));
end
y=[(df-s)';-h];
