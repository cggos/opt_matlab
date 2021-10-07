function afish=fishfollow(fun,afish0,afish1,LB,UB,best_x,parameter)   %×·Î²º¯Êý
[fishnum,m]=size(afish1);
n=0;
f_max=-inf;
max_i=1;
for i=1:fishnum
    if ~isequal(afish1(i,:),afish0)
      if (fishdstc(afish1(i,:),afish0)<parameter.visual)
        n=n+1;
        if fun(afish1(i,:))>f_max
            f_max=fun(afish1(i,:));
            max_i=i;
        end
      end
    end
end
if ((f_max/n)>(parameter.delta*fun(afish0)))
    %r_step=abs(1-(fun(afish0)/f_max))*parameter.step;
    r_step=parameter.step;
    afish=afish0+r_step.*rand(1,m).*(afish1(max_i,:)-afish0+best_x-afish0)/norm(afish1(max_i,:)-2.*afish0+best_x);
    afish=boundtest(afish,LB,UB);
    return
end
afish=fishprey(fun,afish0,LB,UB,best_x,parameter);