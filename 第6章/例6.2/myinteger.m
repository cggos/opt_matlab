function y=myinteger(x)  %�ж�һ����ֵ�Ƿ�������
n=length(x);
y1=ones(1,n);
for i=1:n
   if abs(round(x(i))-x(i))>1.0e-7
       y1(i)=0;
   end
end
a=find(y1==0);
if isempty(a)
    y=[];
else
    y=a;
end