function [G,id1]=graphpts(M)       %ͼ���ж�������Ӿ���ת�Ƴɱߵ����Ӿ���
A=triu(M);         %�����Ӿ����������
num=length(find(A~=0));    %�ߵ���Ŀ
G=zeros(num,num);         %��ת���ɵ������Ӿ���
[id(:,1),id(:,2)]=find(A~=0);
id(:,3)=(1:num)';
[id(:,1),b]=sort(id(:,1));
id(:,2)=id(b,2);
side=cell(1,num);
for i=1:num
    B1=redu(id,i,'r');   %��i����
    [a,c]=find(B1(:,1:2)==id(i,1));
    [b,d]=find(B1(:,1:2)==id(i,2)); 
    side{i}=[B1(a,3);B1(b,3)];
    G(i,side{i})=1;
end
id1=id(:,1:2);