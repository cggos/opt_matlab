function y=myclass(x,m_center)   %根据聚类中心，对样品进行分类
n=size(x,1);
for i=1:n
    d=fy(x(i,:),m_center);
    [val,y(i)]=min(d);
end
