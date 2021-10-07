function afish=fishevaluate(fun,afish0,afish1,LB,UB,best_x,parameter)   %ÓãµÄÆÀ¼Û
afish1=fishfollow(fun,afish0,afish1,LB,UB,best_x,parameter);
afish2=fishswarm(fun,afish0,afish1,LB,UB,best_x,parameter);
afish3=fishprey(fun,afish0,LB,UB,best_x,parameter);
af_best=afish1;
if fun(afish2)>fun(af_best)
    af_best=afish2;
end
if fun(afish3)>fun(af_best)
    af_best=afish3;
end
if fun(af_best)>fun(afish0)
    afish=af_best;
else
   afish=fishmove(afish0,parameter,LB,UB);
end