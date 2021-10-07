function d=fy(x,m_center)  %求样本到类中心的距离
[num,c]=size(m_center);
d=zeros(1,num);
for j=1:num
    for k=1:c
       d(j)=d(j)+(x(k)-m_center(j,k))^2;  %各个样本到各类的距离
    end
end