function afish=fishprey3_1(fun,afish0,parameter)
NC=length(afish0);
for i=1:parameter.try_number
    afish1=afish0;
    afish=afish0;
    n=randperm(NC);
    n=n(1:parameter.visual);
    afish1(n)=abs(afish0(n)-1);
    m=randperm(parameter.visual);
    if length(m)>parameter.step
        m=m(1:parameter.step);
    end
    y1=fun(afish1);
    y2=fun(afish0);
    if y1>y2
       afish(n(m))=afish1(n(m));
       return
    end
end
afish=fishmove3_1(afish0,parameter);