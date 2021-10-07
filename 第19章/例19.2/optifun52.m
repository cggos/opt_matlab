function y=optifun52(x)
NC=length(x);
y1=0;
y=1;
for i=1:NC
    y1=y1+abs(x(i));
    y=y*abs(x(i));
end
y=-y-y1;