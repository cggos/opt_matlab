function y=optifun96(x)
NC=length(x);
y1=0;y2=0;
for i=1:NC
    y1=y1+x(i)^2;
    y2=y2+0.5*i*x(i);
end
y=y1+y2^2+y2^4;