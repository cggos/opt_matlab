function y=optifun73(x)
%y=sin(sqrt(x(1)^2+x(2)^2))-cos(sqrt(abs(x(1)^2-x(2)^2)))-0.02*(x(2)-4.96)^2-0.02*(x(1)-5.87)^2;
y=exp(-(x(1)-0.1)^2)*(sin(5*pi*x(1)^(3/4)))^6;