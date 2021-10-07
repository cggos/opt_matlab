function y=optifun91(x)
NC=length(x);
y=0;
for i=1:NC
    y=y+1/(i+(x(i)-1)^2);
end
y=1/(0.01+y);
