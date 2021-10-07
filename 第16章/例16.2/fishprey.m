function afish=fishprey(fun,afish0,LB,UB,best_x,parameter)    %цыйЁ
n=length(afish0);
for i=1:parameter.try_number
    afish_next=afish0+parameter.visual.*rand(1,n);
    yj=fun(afish_next);
    yi=fun(afish0);
    if yi<yj
      %r_step=abs(1-yi/yj)*parameter.step;
       r_step=parameter.step;
       afish=afish0+r_step.*rand(1,n).*(afish_next-afish0+best_x-afish0)/norm(afish_next-2.*afish0+best_x);
       afish=boundtest(afish,LB,UB);
       return
    end
end
afish=fishmove(afish0,parameter,LB,UB);
