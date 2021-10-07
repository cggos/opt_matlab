function y=progturn(A,b,sigma,f0,k,l)   %线性规划的转轴函数
[r,c1]=size(A);
zuyan=A(k,l);
for i=1:r
   for j=1:c1 
       if i==k
           y1(i,j)=A(k,j)/zuyan;
       else           
           y1(i,j)=A(i,j)-A(k,j)*A(i,l)/zuyan;
       end
   end
end
for i=1:r
    if i==k
        y2(i,1)=b(k)/zuyan;
    else           
        y2(i,1)=b(i)-b(k)*A(i,l)/zuyan;
    end 
end
for i=1:c1
    y3(i)=sigma(i)-A(k,i)*sigma(l)/zuyan;
end
if isempty(f0)
    y4=nan;
else
   y4=f0-b(k)*sigma(l)/zuyan;
end
y=[y1 y2;y3 y4];
   