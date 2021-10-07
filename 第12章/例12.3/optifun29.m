function y=optifun29(x)
data_x=[9 14 21 28 42 57 63 70 79];
data_y=[8.93 10.80 18.59 22.33 39.35 56.11 61.73 64.62 67.08];
n=length(data_x);
y=0;
for i=1:n
   %y=y+(data_y(i)-x(1)*exp(-exp(x(2)-x(3)*data_x(i))))^2;   %Gompertz模型
   y=y+(data_y(i)-x(1)/(1+exp(x(2)-x(3)*data_x(i))))^2;      %Logistic模型
   %y=y+(data_y(i)-x(1)/(1+exp(x(2)-x(3)*data_x(i)))^(1/x(3)))^2;　　%Richards模型
   %y=y+(data_y(i)-x(1)+x(2)*exp(-x(3)*data_x(i)^x(4)))^2;　　　%Weibull模型
end