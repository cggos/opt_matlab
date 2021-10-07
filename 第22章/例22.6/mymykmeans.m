function [label,zhongxin]=mymykmeans(data,K,nc,iter_max)
%nc开始的聚类中心
[m,n]=size(data);%求data的大小，m代表样本的个数，n代表数据特征数
if nargin==2
   iter_max=500;   %设置最大迭代次数
   nc=[];
   zhongxin=zeros(K,n);%存储中心点
end
if nargin==3
    iter_max=500;
    zhongxin=nc;
end
juli=zeros(m,K);%存储各个样本到中心点的距离，行是所有的样本，列是每个样本与中心点的距离
juli1=zeros(m,K);
label=zeros(m,1);%存储各个样本的标签
if isempty(nc)
   for i=1:n
      ma(i)=max(data(:,i));    %每一维最大的数每一列的最大值
      mi(i)=min(data(:,i));    %每一维最小的数
      for j=1:K
          zhongxin(j,i)=ma(i)+(mi(i)-ma(i))*rand();  %随机初始化，不过还是在每一维[min max]中初始化好些
      end      
   end
   %for i=1:K
   %   zhongxin(i,:)=data(ceil(rand*m),:);%随机产生K个聚类中心
   %end
end
 

%%下面开始进行迭代
for diedai=1:iter_max
   for i=1:m
       for j=1:K
          juli1(i,j)=sqrt(sum((data(i,:)-zhongxin(j,:)).^2));%一个样本与所有聚类中心的距离
       end    %一行代表一个样本，K列代表与K个聚类中心的距离
   end
   for i=1:m
      [julisort,zuobiao]=sort(juli1(i,:));%将距离按照由小到大排序
      label(i,1)=zuobiao(1,1);
      %for k=1:K
      %    if zuobiao(1,1)==k%如果最小的距离的标签是k
      %       label(i,1)=k;%则第k个样本的标签就是k
      %    end
      %end
   end
   sumaver=zeros(K,n);%初始化中心点的特征向量sumaver
   geshu=zeros(K,1);%为求中心点的特征向量的平均值而设定的计数器
   %计算聚类中心
   for i=1:m
       sumaver(label(i,1),:)=sumaver(label(i,1),:)+data(i,:);
       geshu(label(i,1),1)=geshu(label(i,1),1)+1;             
   end
   for k=1:K
       sumaver(k,:)=sumaver(k,:)./geshu(k,1);
       zhongxin(k,:)=sumaver(k,:);%更新聚类中心
   end
   myerror=0;
   for i=1:m
      for k=1:K%如果距离不再变化，即中心点不再变化，或者达到了最大的迭代次数，则停止迭代
         myerror=myerror+sum((juli1(i,k)-juli(i,k)).^2);
      end
   end
   juli=juli1;
   juli1=zeros(m,K);
   if myerror==0%如果所有的中心点不再移动
      break;
   end
end
label=label';
