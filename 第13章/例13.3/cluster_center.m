function m=cluster_center(x,pattern) %���������
%patternΪģʽ��Ϊһ���У���ʾÿ����Ʒ��Ӧ�����ţ���[1 2 3 1 2 3 ...]��centernumΪ�����
[n,c]=size(x);
class2=unique(pattern);
centernum=length(class2);
m=zeros(centernum,c);
if n==1
    m(class2,:)=x(1,:);
end
for j=1:centernum
  a=find(pattern==class2(j)); 
  if length(a)==1
       m(j,:)=x(a,:);
  else
       m(j,:)=mean(x(find(pattern==class2(j)),:));
  end
end