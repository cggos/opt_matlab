function [path,dist,dist1]=Mindist(M,varargin)    %ͼ���������·
    %dist      ����̾������
    %dist1     �����ж����ŵ���̾������
    %M             ��ͼ��Ȩֵ����
    %k1,k2�����������������յ�
    %j1,j2         :�ؾ����ĵ�
%��1������ͼ��Ȩֵ�������������������̾���
if nargin==1
    [path,dist,dist0]=pointdist(M);  
elseif nargin==4 && varargin{3}==1     %��k1��k2֮������·
     k1=varargin{1};k2=varargin{2};
     [path,dist,dist0]=pointdist(M,k1,k2); 
elseif nargin==4 && varargin{3}==2     %��k1��k2֮��Ĵζ�·
     k1=varargin{1};k2=varargin{2};
     [path1,distance1,dist0]=pointdist(M,k1,k2);
     n=length(path1); %k1��k2֮������·�ϵĶ�����
     distance2=inf; %����distance2��ֵΪ�����
     for i=1:n-1
         M1=M;
         %ɾ�����·�ϵ�һ����
         M1(path1(i),path1(i+1))=inf;
         M1(path1(i+1),path1(i))=inf;
         [path,distance]=pointdist(M1,k1,k2); %����ͼ�����·���������·
         if distance<distance2
            distance2=distance;
            path2=path;
         end
     end
     path=[{path1} {path2}];
     dist=[{distance1} {distance2}];
elseif nargin==5                                 %����������յ㼰�ؾ��ĵ�
    k1=varargin{1};k2=varargin{2};
    j1=varargin{3};j2=varargin{4};
    [path1,distance1,dist0]=pointdist(M,k1,j1);
    [path2,distance2]=pointdist(M,j1,j2); %����j1-j2�����·��
    [path3,distance3]=pointdist(M,j2,k2); %����j2-k2�����·��
    path1_distance=distance1+distance2+distance3;
    [path4,distance4]=pointdist(M,k1,j2); %����k1-j2�����·��
    [path5,distance5]=pointdist(M,j2,j1); %����j2-j1�����·��
    [path6,distance6]=pointdist(M,j1,k2); %����j1-k2�����·��
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
dist1=zeros(n+1,n+1);   %�������ͼ�Ķ����Ŵ����̾������ľ���
for i=1:n+1
    dist1(1,i)=i-1;     %��һ�дӵڶ���Ԫ�ؿ�ʼ����Ϊͼ�Ķ�����
    dist1(i,1)=i-1;     %��һ�дӵڶ���Ԫ�ؿ�ʼ����Ϊͼ�Ķ�����
end
dist1(1,1)=NaN;         %�������еĵ�һ��Ԫ����NaN�滻
for i=2:n+1
   for j=2:n+1
       dist1(i,j)=dist0(i-1,j-1);
   end
end


