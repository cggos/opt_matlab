function y=boundtest(x,LB,UB,type)    %�߽���
c=size(LB,2);     %���ǵ����������������ļ��,��c��
%ÿһά�����߽���
if nargin==3
    type=1;
end
[m,n]=size(x);   %mΪ��Ʒ��
y=zeros(m,n);
if c==1  %ֻ��һ��
     for i=1:m
       for j=1:n
          if type==1
             p=x(i,j)-UB(j,1);
             q=x(i,j)-LB(j,1);
             if x(i,j)>UB(j,1)
                 y(i,j)=LB(j,1)+p*(UB(j,1)-LB(j,1))/q;
             elseif x(i,j)<LB(j,1)
                 y(i,j)=LB(j,1)+q*(UB(j,1)-LB(j,1))/p;
             else
                y(i,j)=x(i,j);
             end
          elseif type==2
             if x(i,j)<LB(j,1)
                 y(i,j)=min(UB(j,1),2*LB(j,1)-x(i,j));
             elseif x(i,j)>UB(j,1)
                 y(i,j)=max(2*UB(j,1)-x(i,j),LB(j,1));
             else
                 y(i,j)=x(i,j);
             end
          end
       end
    end
else      %��ʱΪ�ṹ����ʽ 
    for i=1:m     %mΪά��
        for j=1:n    %Ϊ�����
            if type==1
               p=x(i,j)-UB(i,j);
               q=x(i,j)-LB(i,j);
               if x(i,j)>UB(i,j)
                   y(i,j)=LB(i,j)+p*(UB(i,j)-LB(i,j))/q;
               elseif x(i,j)<LB(i,j)
                   y(i,j)=LB(i,j)+q*(UB(i,j)-LB(i,j))/p;
               else
                   y(i,j)=x(i,j);
               end
            elseif type==2
               if x(i,j)<LB(i,j)
                 y(i,j)=min(UB(i,j),2*LB(i,j)-x(i,j));
               elseif x(i,j)>UB(j,1)
                 y(i,j)=max(2*UB(i,j)-x(i,j),LB(i,j));
               else
                 y(i,j)=x(i,j);
               end
            end
       end
    end
end
