function y=optifun85(x)
NC=length(x);
y=0;
for i=1:NC
    y=y+(x(i)+0.5)^2;
end
