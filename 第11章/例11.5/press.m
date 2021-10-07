function [y1,b]=press(x,y)
x1=[ones(length(y),1) x];
[b,bint,r]=regress(y,x1,0.01);
y1=0;
for i=1:length(y)
    hi(i)=x1(i,:)*inv(x1'*x1)*x1(i,:)';
    y1=y1+(r(i)/(1-hi(i)))^2;
end
