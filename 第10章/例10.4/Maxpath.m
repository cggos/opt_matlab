function varargout=Maxpath(M,varargin)
    %M                ：网络图的完好概率权值矩阵
    %k1               ：需要求最大可靠性路径的起点
    %k2               ：需要求最大可靠性路径的终点
    %path             ：所给定两点间的最大可靠性路径
    %probability      ：最大可靠性路径的概率
    %capacity         :最大可靠性路径的容量
    %flag             ：指标参数，flag=1表示未找到解
if nargin==3
    k1=varargin{1};k2=varargin{2};
    [m,n]=size(M);
    flag=0;
    M1=zeros(m,n);
    %对完好概率矩阵进行转换
    for i=1:m
      for j=1:n
        if M(i,j)>0 && M(i,j)<1
            M1(i,j)=-log(M(i,j));
        elseif M(i,j)==0
            M1(i,j)=inf;
        end
      end
    end
    %调用最短路算法程序
    [path,distance]=pointdist(M1,k1,k2);  
    if distance<inf
       probability=1;
       for i=1:(length(path)-1)
           probability=probability*M(path(i),path(i+1));  %计算最大可靠性路的概率
       end
    else
       path=0;
       probability=0;
       flag=1;
    end
    %如果没有解，输出如下提示性语句
    if flag==1
       fprintf('顶点%i与顶点%i之间无最大可靠性路：',k1,k2)
    end
    varargout={path,probability,flag};
else
   C=varargin{1};k1=varargin{2};k2=varargin{3};
   flag=0;                      %计算最大可靠性路的默认选择
   capacity = 0;                %给期望最大可靠性容量路的容量赋初值
   k=1;
   %第一步
   while flag==0 && k<100
       [path1,probability1,flag]=Maxpath(M,k1,k2); %调用求最大可靠性路的程序
       if flag==0
        %计算最大可靠性路的容量，以瓶颈弧上的容量为准
          c1=inf;
          for i=1:(length(path1)-1)
               c2=C(path1(i),path1(i+1));
               if c1>c2
                   c1=c2;
               end
          end
          capacity1=c1*probability1; %计算路的期望容量
        %第二步
          C(C<c1)=0;
          M(C<c1)=0;
        %第三步
          if  capacity1>capacity
              capacity=capacity1;
              path=path1;
              probability = probability1;
          end
          k=k+1;
       end
   end
   %如果没有解，输出如下提示性语句
   varargout={path,probability,capacity,flag};
   if flag==1
      fprintf('顶点%i与顶点%i之间无期望最大可靠性容量路：',k1,k2)
   end
end