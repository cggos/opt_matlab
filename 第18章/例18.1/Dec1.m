function pop=Dec1(pop,LB,UB,CodeLen,n)           %二进制转换为十进制编码,n为变量域的序号
NC=size(LB,1);
popsize=size(pop,2);
names=fieldnames(pop);
for i=1:popsize
    pop(i).x=getfield(pop(i),names{n});
end
sublen=CodeLen/NC;
pow_two=2.^(0:sublen)';
maxintval=((2^sublen))-1;
range=UB'-LB';
for i=1:popsize
   start=1;
   fin=sublen;
   for j=1:NC
      tvars(1:sublen)=pop(i).x(start:fin);
      start=start+sublen;
      fin=fin+sublen;
      temp1=0;
      for k=1:sublen
         temp1=temp1+pow_two(k)*tvars(sublen-k+1);
      end
      pop(i).var(j)=(range(j)*(temp1/maxintval))+LB(j);
   end
end
