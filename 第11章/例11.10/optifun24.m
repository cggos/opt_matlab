function y=optifun24(x)
y1=sin(x(1)+x(2))-6*exp(x(1))*x(2);
y2=5*x(1)^2-4*x(2)-100;
y=y1^2+y2^2;
y=-1/(1+sqrt(y));
