function y=isin(x,A)   %ÅĞ¶ÏÊÇ·ñÔÚAÖĞ
k=0;
for i=1:length(A)
   if abs(x-A(i))<=0.0001
      k=k+1;
      y(k)=i;
      break
   else
       y=0;
   end
end
