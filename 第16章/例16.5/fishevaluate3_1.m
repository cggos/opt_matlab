function afish=fishevaluate3_1(fun,afish0,afish1,parameter)
af_follow=fishfollow3_1(fun,afish0,afish1,parameter);
af_swarm=fishswarm3_1(fun,afish0,afish1,parameter);
af_prey=fishprey3_1(fun,afish0,parameter);
afish2=af_follow;
if fun(af_swarm)<fun(afish2)
    afish2=af_swarm;
end
if fun(af_prey)<fun(afish2)
    afish2=af_prey;
end
if fun(afish2)<fun(afish0)
    afish=afish2;
else
    afish=fishmove3_1(afish0,parameter);
end