function y=dis2(x,m,pattern) %�����ھ���,mΪ�������ģ�patternΪ���ģʽ
class=size(m,1);    %���
for j=1:class
    d(j)=0;
    for i=1:size(x,1)
       if pattern(i)==j
           d(j)=d(j)+norm(x(i,:)-m(j,:));
       end
    end
end
y=sum(d);
