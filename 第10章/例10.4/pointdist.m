function [path,distance,dist0]=pointdist(M,k1,k2)   %图中任意两点间的最短路
n=length(M);      %图的定点数
if nargin==1
    k1=1;
    k2=n;
end
dist0=M;          %给距离矩阵赋初值
m=1;
while m<=n
    for i=1:n
        for j=1:n
            if dist0(i,j)>dist0(i,m)+dist0(m,j)
                dist0(i,j)=dist0(i,m)+dist0(m,j);
            end
        end
    end
    m=m+1;
end
distance=dist0(k1,k2);
Path0=zeros(1,n);
k=1;
Path0(k)=k2;
kk=k2;
V=inf*ones(1,n);
while kk~=k1
    for i=1:n
       V(1,i)=dist0(k1,kk)-M(i,kk);
       if V(1,i)==dist0(k1,i)
            Path0(k+1)=i; kk=i; k=k+1;
       end
    end
end
k=1;
wrow=find(Path0~=0);
for j=length(wrow):-1:1
    path(k)=Path0(wrow(j));
    k=k+1;
end
