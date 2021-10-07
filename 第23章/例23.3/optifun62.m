function y=optifun62(x)
a=0;
for i=1:3
    a=a+abs(x(i))^0.8+5*sin(x(i)^3);
end
b=0;
for i=1:2
  b=b+(-10*exp(-0.2*sqrt(x(i)^2+x(i+1)^2)));
end
y=(x(4)*a+x(5)*b+5000*(x(4)+x(5)-1)^2);
y=-y;
