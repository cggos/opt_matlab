function [y1,y2]=optifun26(x)
x=[0 x 1];
c=length(x);
for j=1:c
   y3(j)=my_fun(x(j));
end
for j=1:c-1
     d(j)=x(j+1)-x(j);
     dmid(j)=(x(j+1)+x(j))/2;
     y_mid(j)=my_fun(dmid(j));
     y_max(j)=max([y3(j),y3(j+1),y_mid(j)]);
     y_min(j)=min([y3(j),y3(j+1),y_mid(j)]);
     dd1(j)=(y_max(j)-y_min(j))^2*d(j);
     dd2(j)=y_mid(j)*d(j);
end
y1=sum(dd1)/2;%适应度函数值
y2=sum(dd2);  %定积分值



function y=my_fun(x)
y=x*sin(x)*sin(100*x);
