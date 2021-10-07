function y=optifun40(x,data,target)
a=findzeros(x);
if strcmp(a,'all')
    y=inf;
else
  N=size(data,2);
  data=data(:,find(x==1));
  m=cluster_center(data,target);%求各类中心
  y1=cluster_dis(data,target,m);    %类间、类内距离
  y2=crossvalindtest(data,10,3,target);   %正确率
  %knnmd1=fitcknn(data,target','NumNeighbors',5,'Standardize',1);
  %CVKNNMdl=crossval(knnmd1);
  %y2=1-kfoldLoss(CVKNNMdl);
  y=y1-(y2-0.1*length(find(x==1))/N);
end
