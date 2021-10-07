function y=findcount(x)      %从一组数中找每个数出现的次数
num=length(x);
x1=unique(x);
n2=length(x1);
y=zeros(n2,2);
for i=1:n2
   y(i,1)=x1(i);
   for j=1:num
       if x1(i)==x(j)
           y(i,2)=y(i,2)+1;
       end
   end
end
[a,b]=sort(y(:,2));
y(:,1)=y(b,1);
y(:,2)=a;