function [label,zhongxin]=mymykmeans(data,K,nc,iter_max)
%nc��ʼ�ľ�������
[m,n]=size(data);%��data�Ĵ�С��m���������ĸ�����n��������������
if nargin==2
   iter_max=500;   %��������������
   nc=[];
   zhongxin=zeros(K,n);%�洢���ĵ�
end
if nargin==3
    iter_max=500;
    zhongxin=nc;
end
juli=zeros(m,K);%�洢�������������ĵ�ľ��룬�������е�����������ÿ�����������ĵ�ľ���
juli1=zeros(m,K);
label=zeros(m,1);%�洢���������ı�ǩ
if isempty(nc)
   for i=1:n
      ma(i)=max(data(:,i));    %ÿһά������ÿһ�е����ֵ
      mi(i)=min(data(:,i));    %ÿһά��С����
      for j=1:K
          zhongxin(j,i)=ma(i)+(mi(i)-ma(i))*rand();  %�����ʼ��������������ÿһά[min max]�г�ʼ����Щ
      end      
   end
   %for i=1:K
   %   zhongxin(i,:)=data(ceil(rand*m),:);%�������K����������
   %end
end
 

%%���濪ʼ���е���
for diedai=1:iter_max
   for i=1:m
       for j=1:K
          juli1(i,j)=sqrt(sum((data(i,:)-zhongxin(j,:)).^2));%һ�����������о������ĵľ���
       end    %һ�д���һ��������K�д�����K���������ĵľ���
   end
   for i=1:m
      [julisort,zuobiao]=sort(juli1(i,:));%�����밴����С��������
      label(i,1)=zuobiao(1,1);
      %for k=1:K
      %    if zuobiao(1,1)==k%�����С�ľ���ı�ǩ��k
      %       label(i,1)=k;%���k�������ı�ǩ����k
      %    end
      %end
   end
   sumaver=zeros(K,n);%��ʼ�����ĵ����������sumaver
   geshu=zeros(K,1);%Ϊ�����ĵ������������ƽ��ֵ���趨�ļ�����
   %�����������
   for i=1:m
       sumaver(label(i,1),:)=sumaver(label(i,1),:)+data(i,:);
       geshu(label(i,1),1)=geshu(label(i,1),1)+1;             
   end
   for k=1:K
       sumaver(k,:)=sumaver(k,:)./geshu(k,1);
       zhongxin(k,:)=sumaver(k,:);%���¾�������
   end
   myerror=0;
   for i=1:m
      for k=1:K%������벻�ٱ仯�������ĵ㲻�ٱ仯�����ߴﵽ�����ĵ�����������ֹͣ����
         myerror=myerror+sum((juli1(i,k)-juli(i,k)).^2);
      end
   end
   juli=juli1;
   juli1=zeros(m,K);
   if myerror==0%������е����ĵ㲻���ƶ�
      break;
   end
end
label=label';
