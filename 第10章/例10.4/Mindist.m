function [path,dist,dist1]=Mindist(M,varargin)    %图论中求最短路
    %dist      ：最短距离矩阵
    %dist1     ：带有顶点编号的最短距离矩阵
    %M             ：图的权值矩阵
    %k1,k2　　　　　：起点和终点
    %j1,j2         :必经过的点
%（1）根据图的权值矩阵求解任意两点间的最短距离
if nargin==1
    [path,dist,dist0]=pointdist(M);  
elseif nargin==4 && varargin{3}==1     %求k1、k2之间的最短路
     k1=varargin{1};k2=varargin{2};
     [path,dist,dist0]=pointdist(M,k1,k2); 
elseif nargin==4 && varargin{3}==2     %求k1、k2之间的次短路
     k1=varargin{1};k2=varargin{2};
     [path1,distance1,dist0]=pointdist(M,k1,k2);
     n=length(path1); %k1、k2之间的最短路上的顶点数
     distance2=inf; %赋予distance2初值为无穷大
     for i=1:n-1
         M1=M;
         %删除最短路上的一条边
         M1(path1(i),path1(i+1))=inf;
         M1(path1(i+1),path1(i))=inf;
         [path,distance]=pointdist(M1,k1,k2); %求新图的最短路——次最短路
         if distance<distance2
            distance2=distance;
            path2=path;
         end
     end
     path=[{path1} {path2}];
     dist=[{distance1} {distance2}];
elseif nargin==5                                 %经过起点与终点及必经的点
    k1=varargin{1};k2=varargin{2};
    j1=varargin{3};j2=varargin{4};
    [path1,distance1,dist0]=pointdist(M,k1,j1);
    [path2,distance2]=pointdist(M,j1,j2); %计算j1-j2的最短路径
    [path3,distance3]=pointdist(M,j2,k2); %计算j2-k2的最短路径
    path1_distance=distance1+distance2+distance3;
    [path4,distance4]=pointdist(M,k1,j2); %计算k1-j2的最短路径
    [path5,distance5]=pointdist(M,j2,j1); %计算j2-j1的最短路径
    [path6,distance6]=pointdist(M,j1,k2); %计算j1-k2的最短路径
    path2_distance=distance4+distance5+distance6;
    if path1_distance<path2_distance
       dist=path1_distance; 
       path=[path1 path2(2:length(path2)) path3(2:length(path3))];
    else
       dist=path2_distance;
       path=[path4 path5(2:length(path5)) path6(2:length(path6))];
    end
end
n=length(M);
dist1=zeros(n+1,n+1);   %定义带有图的顶点编号存放最短距离结果的矩阵
for i=1:n+1
    dist1(1,i)=i-1;     %第一行从第二个元素开始依次为图的顶点编号
    dist1(i,1)=i-1;     %第一列从第二个元素开始依次为图的顶点编号
end
dist1(1,1)=NaN;         %将矩阵中的第一个元素以NaN替换
for i=2:n+1
   for j=2:n+1
       dist1(i,j)=dist0(i-1,j-1);
   end
end


