function y=optifun45(x)
I=[0.10 0.30 0.50 1
    0.17 0.33 0.50 1
    0.32 0.48 0.48 1
    0.12 0.30 0.50 1
    0.10 0.30 0.50 1];
PI0=[0.2564 0.3509 0.4803 0.900];
for i=1:4
    y1=0;
    for j=1:5
        y1=y1+I(j,i);
    end
    y1=y1/5;
    y=(PI0(i)-x(1)*y1^x(2))^2;
end
y=-y;
    
