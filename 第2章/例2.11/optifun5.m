function [f,g]=optifun5(x)
f=3*x(1)^2+2*x(1)*x(2)+x(2)^2;    
if nargout>1
   g(1)=6*x(1)+2*x(2);   %Ìİ¶È
   g(2)=2*x(1)+2*x(2);
end