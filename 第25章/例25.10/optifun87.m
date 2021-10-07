function y=optifun87(x)
%y=sin(sqrt(x(1)^2+x(2)^2))-cos(sqrt(abs(x(1)^2-x(2)^2)))-0.02*(x(2)-4.96)^2-0.02*(x(1)-5.87)^2;
%y=-y;
y=10*sin(5*pi*sqrt(0.01*x))*exp(-(0.01*x-0.5)^2+1)+30;
y=-y;
