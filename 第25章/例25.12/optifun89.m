function y=optifun89(x)
%y=20*exp(0.2*sqrt((x(1)^2+x(2)^2)/2))+exp((cos(2*pi*x(1))+cos(2*pi*x(2)))/2);
%f(0,0)=22.7128
y=(x(1)^2+x(2)^2)/2-cos(2*x(1))*cos(2*x(2));
