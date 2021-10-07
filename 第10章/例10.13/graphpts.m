function [G,id1]=graphpts(M)       %图论中顶点的连接矩阵转制成边的连接矩阵
A=triu(M);         %求连接矩阵的上三角
num=length(find(A~=0));    %边的数目
G=zeros(num,num);         %边转化成点后的连接矩阵
[id(:,1),id(:,2)]=find(A~=0);
id(:,3)=(1:num)';
[id(:,1),b]=sort(id(:,1));
id(:,2)=id(b,2);
side=cell(1,num);
for i=1:num
    B1=redu(id,i,'r');   %第i条边
    [a,c]=find(B1(:,1:2)==id(i,1));
    [b,d]=find(B1(:,1:2)==id(i,2)); 
    side{i}=[B1(a,3);B1(b,3)];
    G(i,side{i})=1;
end
id1=id(:,1:2);