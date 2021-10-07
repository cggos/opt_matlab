function afish=fishfollow3_1(fun,afish0,afish1,parameter)
fishnum=size(afish1,1);
n=0;
f_max=-inf;
max_i=1;
for i=1:fishnum
    if (fishdstc3_1(afish1(i,:),afish0)<parameter.visual)
        n=n+1;
        f=fun(afish1(i,:));
       if f>f_max
           f_max=f;
           max_i=i;
       end
    end
end
if (f_max/n>parameter.delta*fun(afish0))
    afish=afish0;
    m=find((afish1(max_i,:)-afish0)~=0);     
    n=randperm(length(m));
    if length(m)>parameter.step
       n=n(1:parameter.step);
    end
    afish(m(n))=afish1(max_i,m(n));
    return
else
    afish=fishprey3_1(fun,afish0,parameter);
end
