function d=fy(x,m_center)  %�������������ĵľ���
[num,c]=size(m_center);
d=zeros(1,num);
for j=1:num
    for k=1:c
       d(j)=d(j)+(x(k)-m_center(j,k))^2;  %��������������ľ���
    end
end