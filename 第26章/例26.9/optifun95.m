function y=optifun95(x)
NC=length(x);
y=0;
for i=1:NC
    y1=0;
    for j=1:i
        y1=y1+x(j);
    end
    y=y+y1^2;
end
