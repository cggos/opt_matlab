function y=optifun37(x,data,target)
if all(x==0)
   y=inf;
else
  N=size(data,2);
  data1=data(:,find(x==1));
  m=cluster_center(data1,target);%求各类中心
  y1=cluster_dis(data1,target,m);    %类间、类内距离
  y2=crossvalindtest(data1,6,3,target);   %正确率
  y=y1-(y2-0.1*length(find(x==1))/N);
end

