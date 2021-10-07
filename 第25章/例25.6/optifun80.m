function y=optifun80(x)
y=x(1)*sin(sqrt(abs(x(2)+1-x(1))))*cos(sqrt(abs(x(2)+x(1)+1)))...
    +(x(2)+1)*cos(sqrt(abs(x(2)+1-x(1))))*sin(sqrt(abs(x(2)+x(1)+1)));

