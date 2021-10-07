function varargout=Maxpath(M,varargin)
    %M                ������ͼ����ø���Ȩֵ����
    %k1               ����Ҫ�����ɿ���·�������
    %k2               ����Ҫ�����ɿ���·�����յ�
    %path             �����������������ɿ���·��
    %probability      �����ɿ���·���ĸ���
    %capacity         :���ɿ���·��������
    %flag             ��ָ�������flag=1��ʾδ�ҵ���
if nargin==3
    k1=varargin{1};k2=varargin{2};
    [m,n]=size(M);
    flag=0;
    M1=zeros(m,n);
    %����ø��ʾ������ת��
    for i=1:m
      for j=1:n
        if M(i,j)>0 && M(i,j)<1
            M1(i,j)=-log(M(i,j));
        elseif M(i,j)==0
            M1(i,j)=inf;
        end
      end
    end
    %�������·�㷨����
    [path,distance]=pointdist(M1,k1,k2);  
    if distance<inf
       probability=1;
       for i=1:(length(path)-1)
           probability=probability*M(path(i),path(i+1));  %�������ɿ���·�ĸ���
       end
    else
       path=0;
       probability=0;
       flag=1;
    end
    %���û�н⣬���������ʾ�����
    if flag==1
       fprintf('����%i�붥��%i֮�������ɿ���·��',k1,k2)
    end
    varargout={path,probability,flag};
else
   C=varargin{1};k1=varargin{2};k2=varargin{3};
   flag=0;                      %�������ɿ���·��Ĭ��ѡ��
   capacity = 0;                %���������ɿ�������·����������ֵ
   k=1;
   %��һ��
   while flag==0 && k<100
       [path1,probability1,flag]=Maxpath(M,k1,k2); %���������ɿ���·�ĳ���
       if flag==0
        %�������ɿ���·����������ƿ�����ϵ�����Ϊ׼
          c1=inf;
          for i=1:(length(path1)-1)
               c2=C(path1(i),path1(i+1));
               if c1>c2
                   c1=c2;
               end
          end
          capacity1=c1*probability1; %����·����������
        %�ڶ���
          C(C<c1)=0;
          M(C<c1)=0;
        %������
          if  capacity1>capacity
              capacity=capacity1;
              path=path1;
              probability = probability1;
          end
          k=k+1;
       end
   end
   %���û�н⣬���������ʾ�����
   varargout={path,probability,capacity,flag};
   if flag==1
      fprintf('����%i�붥��%i֮�����������ɿ�������·��',k1,k2)
   end
end