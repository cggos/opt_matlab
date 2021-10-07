function y=optifun50(x)
%y=abs(sin(30*x))*(1-x/2);
y=sin(sum(abs(x-5)))/sum(abs(x-5));