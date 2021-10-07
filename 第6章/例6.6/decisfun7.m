function  u=decisfun7(k,x)
q=10*[6 7 12 6];
if q(k)-x<0
    u=0:100;
else
    u=q(k)-x:100;
end
u=u(:);