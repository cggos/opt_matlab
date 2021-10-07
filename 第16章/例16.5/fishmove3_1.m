function afish=fishmove3_1(afish0,parameter)
NC=length(afish0);
afish=afish0;
if rand<0.8
  m=ceil(NC*rand(1,parameter.step));
  afish(m)=abs(afish0(m)-1);
else
   afish=afish0;
end