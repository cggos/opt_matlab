function y=optifun67(x)
a1=x(1)^2+x(2)^2+x(3)^2-3;
a2=x(1)^2+x(2)^2+x(1)*x(2)+x(1)+x(2)-5;
a3=x(1)+x(2)+x(3)-3;
y=sum(abs([a1 a2 a3]));
y=-y;