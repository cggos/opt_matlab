function [y1,y2]=optifun58(x)   %函数值及导数矩阵
y1(1,1)=4*x(1)^2+x(2)^2-4;
y1(2,1)=x(1)+x(2)-sin(x(1)-x(2));
y2=[8*x(1) 2*x(2)
    1-cos(x(1)-x(2)) 1+cos(x(1)+x(2))];
