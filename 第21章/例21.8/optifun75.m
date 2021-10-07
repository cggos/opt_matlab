function [y,y1]=optifun75(x,a,p,c)
y1=x*p';
if x*a'-c<=0
    y=y1;
else
    y=-1;
end

