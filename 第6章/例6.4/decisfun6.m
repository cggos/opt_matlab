function  [u1,u2]=decisfun6(k,x1,x2)
w=[3 4 5];v=[8 6 4];
u1=0:1:min(x1/w(k),x2/v(k));
u2=1;
