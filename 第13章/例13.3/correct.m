function [y,b]=correct(x1,x2)     %�ж�x1,x2�Ƿ���ȣ�����ȷ���Ƕ���
[r1,c1]=size(x1);
[r2,c2]=size(x2);
if r1~=r2
    a=x1-x2';
elseif c1~=c2
    a=x1-x2';
else
   a=x1-x2; 
end
b=find(a~=0);
y=length(find(a==0))/length(x1);

