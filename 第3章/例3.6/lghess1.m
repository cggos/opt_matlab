function y=lghess1(fun,hfun,x0,mu,x_syms)   %拉格朗日函数的Hess阵
Hf=hessian(fun,x_syms);
dHf=eval(subs(Hf,x_syms,x0));
HH=myhessian1(hfun,x_syms);
num=size(hfun,1);
if num==1
    s=mu*eval(subs(HH{1},x_syms,x0));
else
  s=0;
  for i=1:num
     s=s+mu(i)*eval(subs(HH{i},x_syms,x0));
  end
end
y=dHf-s;
