function y=myclass(x,m_center)   %���ݾ������ģ�����Ʒ���з���
n=size(x,1);
for i=1:n
    d=fy(x(i,:),m_center);
    [val,y(i)]=min(d);
end
