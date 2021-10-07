function [y,idx]=trimatrix(x,type)      %�����Ǿ�����Ԫ�ص���С���󣩼��ͣ�xΪ������
[r,c]=size(x);
if nargin==1
    type='min';   %Ĭ������Сֵ
end
if strcmp(type,'min')
   x1=inf(r,c);
   for i=1:r
      for j=i+1:c
         x1(i,j)=x(i,j);
      end
   end
   y=min(min(x1));
   [idx(1),idx(2)]=find(x1==y);
elseif strcmp(type,'max')
   x1=-inf(r,c);
   for i=1:r
      for j=i+1:c
         x1(i,j)=x(i,j);
      end
   end
   y=max(max(x1));
   [idx(1),idx(2)]=find(x1==y);
elseif strcmp(type,'sum')   %���
    for i=1:r-1
        y(i)=sum(x(i,i+1:end));
    end
end
