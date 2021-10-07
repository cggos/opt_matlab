function y=Dec(LB,UB,x,CodeLen)           %二进制转换为十进制编码
NC=size(LB,1);
sublen=CodeLen/NC;
pow_two=2.^(0:sublen)';
maxintval=((2^sublen))-1;
range=UB-LB;
start=1;
fin=sublen;
for j=1:NC
   tvars(1:sublen)=x(start:fin);
   start=start+sublen;
   fin=fin+sublen;
   temp1=0;
   for k=1:sublen
       temp1=temp1+pow_two(k)*tvars(sublen-k+1);
   end
   y(j)=(range(j)*(temp1/maxintval))+LB(j);
end
