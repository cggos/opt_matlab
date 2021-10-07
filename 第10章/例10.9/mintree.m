function [Tree,Weight]=mintree(M,type1,type2)   %ͼ��������С֧����
    %Tree     ����С�������ıߵļ���
    %Weight   ����С��������Ȩ֮��
    %M        ��ͼ������Ȩֵ�������һ�ֱ�ʾ�������þ�����һ��3��n�ľ���nΪͼ�ı�����ÿ�е�3��Ԫ�طֱ��ʾ�ߵ���㡢�յ��Ȩֵ��
    %flag     �������������Ʋ���
%���ͼ��Ȩֵ����Ϊ1ά�����������ת��
if nargin==2
    type2='min';
end
if strcmp(type1,'k')    %Kruskal�㷨
  y=findzeros(diag(M));
  if strcmp(y,'all')
    n=size(M,2);           %��ͼ�Ķ�����
    m=sum(sum(M~=0))/2;    %��ͼ�ı���
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
  %���ͼ��Ȩֵ����Ϊ3ά����ֱ�ӽ������
  n=max(max(M1(1:2,:)));     %��ͼ�Ķ�����
  m=size(M1,2);              %��ͼ�ı���
  if strcmp(type2,'min')
     [B,i]=sortrows(M1',3);     %��Ȩ�ķǼ�˳���������б�
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
elseif strcmp(type1,'p')     %Prim�㷨
   n=length(M);         %��ͼ�Ķ�����
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