function y=optifun78(x)
NC=length(x);
y=0;
for i=1:NC-1
    y=y+100*(x(i+1)-x(i)^2)^2+(1-x(i))^2;
end