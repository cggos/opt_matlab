function [c,ceq]=mycon1(x)
ceq=x(2)*x(3)+x(1)*x(2)+x(1)*x(3)-75;
c=[];