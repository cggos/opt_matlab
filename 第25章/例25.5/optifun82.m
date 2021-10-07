function y=optifun82(x)
NC=length(x);
y1=0;y2=1;
for i=1:NC
    y1=y1+x(i)^2/4000;
    y2=y2*cos(x(i)/sqrt(i));
end
y=y1-y2+1;