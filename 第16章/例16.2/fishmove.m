function afish=fishmove(afish0,parameter,LB,UB)      %���
n=length(afish0);
afish=afish0+parameter.visual.*unifrnd(-1,1,1,n);
afish=boundtest(afish,LB,UB);
