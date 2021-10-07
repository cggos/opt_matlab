function [Road,Distance]=graphTSP(M,type)   %ͼ���е�TSP����
    %Road        ������·��
    %Distance    ������·�ߵĳ���
    %M           ��ͼ��Ȩֵ���󣨵ط���ط�֮��ľ������
if type==1     %����Ȧ��
  n=size(M,2);                        %��ͼ�Ķ����������е���Ŀ��
  Road=[linspace(1,n,n) 1];          %������ʼ·�ߣ�1��n��1
  %Road = [1 2 5 6 4 3 1];
  Road1=Road;
  if n>3                             %ֻ���ڳ��е���Ŀ����3���Ŵ�����������
    for r=4:n+1
        for i=1:r-3
            for j=i+2:r-1
                if (M(Road(i),Road(j))+M(Road(i+1),Road(j+1))<M(Road(i),Road(i+1))+M(Road(j),Road(j+1)))
                    Road1(1:i)=Road(1:i);
                    for k=i+1:j
                        Road1(k)=Road(j+i+1-k);
                    end
                    Road1(j+1:r)=Road(j+1:r);
                end
            end
        end
    end
  elseif n<=3
    if n<=2
        fprintf('������������·��������')
    else
        fprintf('����·������')
    end
  end
  Road=Road1;
  Distance=0;
  for i=1:n
    Distance=Distance+M(Road(i),Road(i+1));
  end
elseif type==2   %ö�ٷ�
  n=size(M,1);
  r=exhaus(n);
  [num,c]=size(r);
  r1=[r ones(num,1)];
  d=zeros(num,1);
  Distance=inf;
  for i=1:num
    for j=1:c
       d(i)=d(i)+M(r1(i,j),r1(i,j+1));
    end
    if d(i)<Distance
        Distance=d(i);
    end
  end
  Road=r1(find(d==Distance),:);
end

    
    
    