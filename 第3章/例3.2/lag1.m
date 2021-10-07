function y=lag1(fun,hfun,gfun,x0,mu,lam,sigma)
f=fun(x0);
if ~isempty(hfun)
   he=hfun(x0);
   l=length(he);
   y=f;s1=0;
   for i=1:l
     y=y-he(i)*mu(i);
     s1=s1+he(i)^2;
   end
   y=y+0.5*sigma*s1;
else
   y=f;
end
if ~isempty(gfun)
    gi=gfun(x0);
    m=length(gi);
    s2=0;
    for i=1:m
      s3=max(0,lam(i)-sigma*gi(i));
      s2=s2+s3^2-lam(i)^2;
    end
    y=y+s2/(2.0*sigma);
else
    y=y;
end