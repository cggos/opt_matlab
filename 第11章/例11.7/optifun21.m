function y=optifun21(x)
y=0;
for i=1:length(x)
    y=y+x(i)^2-10*cos(2*pi*x(i))+10;
end