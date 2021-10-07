function y=optifun88(x)
%y=20*exp(0.2*sqrt((x(1)^2+x(2)^2)/2))+exp((cos(2*pi*x(1))+cos(2*pi*x(2)))/2);
%f(0,0)=22.7128
y=0.4+sinc(4*x)+1.1*sinc(4*x+2)+0.8*sinc(6*x-2)+0.7*sinc(6*x-4);
y=-y;