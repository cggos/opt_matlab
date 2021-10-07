function y=dis2(x,m,pattern) %求类内距离,m为聚类中心，pattern为类别模式
class=size(m,1);    %类别
for j=1:class
    d(j)=0;
    for i=1:size(x,1)
       if pattern(i)==j
           d(j)=d(j)+norm(x(i,:)-m(j,:));
       end
    end
end
y=sum(d);
