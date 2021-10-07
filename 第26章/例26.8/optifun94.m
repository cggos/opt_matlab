function y=optifun94(x)
NC=length(x);
y=0;
for i=1:NC
    y=y+x(i)^4-16*x(i)^2+5*x(i);
end
y=y/NC;