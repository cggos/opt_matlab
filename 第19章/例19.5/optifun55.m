function y=optifun55(x)
x=realbit(x);
y1=3*x(1)*x(2)*x(4)+2*x(2)*x(3)*x(4)+9*x(2)*x(4)+5*x(1)-2*x(2)+7*x(4)-17;
y2=x(1)*x(2)*x(3)*x(4)+2*x(2)*x(3)*x(4)+x(2)*x(4)+6*x(1)-2*x(2)+6*x(4)-6;
y3=x(1)*x(2)*x(3)*x(4)+2*x(2)*x(3)*x(4)+x(2)*x(4)+7*x(1)-2*x(2)+7*x(4)-7;
y=y1+5000*(y2^2+y3^2);
if y>0
    y=1/(1+y);
else
    y=1+abs(y);
end
