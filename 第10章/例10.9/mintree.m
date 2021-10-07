function [Tree,Weight]=mintree(M,type1,type2)   %图论中求最小支撑树
    %Tree     ：最小生成树的边的集合
    %Weight   ：最小生成树的权之和
    %M        ：图（树）权值矩阵的另一种表示方法，该矩阵是一个3×n的矩阵，n为图的边数；每列的3个元素分别表示边的起点、终点和权值。
    %flag     ：变量个数控制参数
%如果图的权值矩阵为1维，则进行如下转换
if nargin==2
    type2='min';
end
if strcmp(type1,'k')    %Kruskal算法
  y=findzeros(diag(M));
  if strcmp(y,'all')
    n=size(M,2);           %求图的顶点数
    m=sum(sum(M~=0))/2;    %求图的边数
    M1=zeros(3,m);
    k=1;
    for i=1:n
        for j=i+1:n
            if M(i,j)~=0
                M1(1,k)=i;
                M1(2,k)=j;
                M1(3,k)=M(i,j);
                k=k+1;
            end
        end
    end
  else
    M1=M;
  end
  %如果图的权值矩阵为3维，则直接进行求解
  n=max(max(M1(1:2,:)));     %求图的顶点数
  m=size(M1,2);              %求图的边数
  if strcmp(type2,'min')
     [B,i]=sortrows(M1',3);     %按权的非减顺序重新排列边
  elseif strcmp(type2,'max');
      B=-sortrows(-M1',3);
  end
  B=B';
  Tree=[];
  Weight=0;
  k=1;
  q=1:n;
  for i=1:m
    if q(B(1,i))~=q(B(2,i));
        Tree(1:2,k)=B(1:2,i);
        Weight=Weight+B(3,i);
        k=k+1;
        qmin=min(q(B(1,i)),q(B(2,i)));
        qmax=max(q(B(1,i)),q(B(2,i)));
        for j=1:n
            if q(j)==qmax
                q(j)=qmin;
            end
        end
    end
    if k==n
        break;
    end
  end
elseif strcmp(type1,'p')     %Prim算法
   n=length(M);         %求图的顶点数
   M(M==0)=inf;
   k=1:n;
   listV(k)=0;
   listV(1)=1;
   e=1;
   while (e<n)
      if strcmp(type2,'min')
         min1=inf;
      elseif strcmp(type,'max')
          max1=inf;
      end
      for i=1:n
        if listV(i)==1
            for j=1:n
                if strcmp(type2,'min')
                  if listV(j)==0 && min1>M(i,j)
                    min1=M(i,j);
                    B=M(i,j);
                    S=i;
                    D=j;
                  end
                elseif strcmp(type2,'max')
                   if listV(j)==0 && max1>M(i,j)
                    max1=M(i,j);
                    B=M(i,j);
                    S=i;
                    D=j;
                   end
                end
            end
        end
      end
      listV(D)=1;
      distance(e)=B;
      source(e)=S;
      destination(e)=D;
      e=e+1;
   end
   Tree=[source;destination];
   for g=1:e-1
      Weight(g)=M(Tree(1,g),Tree(2,g));
   end
end