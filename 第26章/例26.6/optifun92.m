function y=optifun92(x)
%y=4+4.5*x(1)-4*x(2)+x(1)^2+2*x(2)^2-2*x(1)*x(2)+x(1)^4-2*x(1)^2*x(2);
NC=length(x);
y1=0;y2=1;
for i=1:NC
    y1=y1+abs(x(i));
    y2=y2*abs(x(i));
end
y=y1+y2;
