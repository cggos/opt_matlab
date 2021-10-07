function xmin=mygoldcut1(fun,x0,d0,a,b,x_syms,esp)   %黄金分割法求最小点
if nargin==6
    esp=1e-6;
end
n=0;
while (abs((a-b))>esp)
    if n>500
       error('无解');
    end
    n=n+1;
    x1=a+(b-a)*0.618;
    if isa(fun,'function_handle')
        y1=fun(x0+x1*d0');
    else
        y1=eval(subs(fun,x_syms,x0+x1*d0'));
    end
    x2=b-(b-a)*0.618;
    if isa(fun,'function_handle')
        y2=fun(x0+x2*d0');
    else
        y2=eval(subs(fun,x_syms,x0+x2*d0'));
    end
    if y1<y2
        a=x2;
    elseif y1>y2
        b=x1;
    elseif y1==y2
        a=x2;
        b=x1;
    end  
end
xmin=(a+b)/2;