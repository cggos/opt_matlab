function y=matrixinsert(x,x1,b)   %����x1��ָ���У��У�,x1��Ԫ��
if iscell(x1)
    n=length(x1);
   [b,dex]=sort(b,'descend');   %�Ӵ�С����
   x1=x1(dex);
   for i=1:n
      y=[x(1:b(i)) x1{i} x(b(i)+1:end)];
      x=y;
   end
else
    n=length(x1);
    [b,dex]=sort(b,'descend');   %�Ӵ�С����
    x1=x1(dex);
    for i=1:n
       y=[x(1:b(i)) x1(i) x(b(i)+1:end)];
       x=y;
    end
end
    