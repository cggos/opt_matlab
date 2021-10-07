function y=matrixinsert(x,x1,b)   %插入x1到指定行（列）,x1是元胞
if iscell(x1)
    n=length(x1);
   [b,dex]=sort(b,'descend');   %从大到小排列
   x1=x1(dex);
   for i=1:n
      y=[x(1:b(i)) x1{i} x(b(i)+1:end)];
      x=y;
   end
else
    n=length(x1);
    [b,dex]=sort(b,'descend');   %从大到小排列
    x1=x1(dex);
    for i=1:n
       y=[x(1:b(i)) x1(i) x(b(i)+1:end)];
       x=y;
    end
end
    