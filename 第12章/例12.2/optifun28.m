function y=optifun28(x)
a=0;b=0;
for i=1:5
    a=a+i*cos((i+1)*x(1)+i);
    b=b+i*cos((i+1)*x(2)+i);
end
y=a*b;
