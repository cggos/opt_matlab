function y=optifun22(x)
c=[160 87 18 71 176 101  35 145 117 54];
a=[198 30 167 130 35 20 105 196 94 126];
b=546;
a1=find(x==1);
if sum(a(a1))>b
   y=-1;
else
   y=sum(c(a1));
end
