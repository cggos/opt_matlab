function [y,y2]=exhaus(x)        %��ٷ��о����еĿ���
if iscell(x)&&~ischar(x{1})
  r=size(x,1);       %xһ��ΪԪ��
  a=1;
  for i=1:r
    n(i)=length(x{i});
    a=a*n(i);     %���е�ȫ�������
    num(i)=a;
  end
  y=zeros(a,r);
  y(1:n(1),1)=x{1}(1:end)';
  for i=2:r
    y(1:num(i),1:i-1)=repmat(y(1:num(i-1),1:i-1),n(i),1);
    for j=1:n(i)
        y((j-1)*num(i-1)+1:j*num(i-1),i)=x{i}(j)*ones(num(i-1),1);
    end
  end
elseif length(x)==1
    x=1:x;
    n=length(x);   %x���е���������
    y=ones(factorial(n-1),n);
    y(1:n-1,1)=x(1);
    y(1:n-1,2)=x(2:end);
    y1=y(1:n-1,1:2);
    for i=3:n
        m=n-(i-1);    %ʣ�������    
        n1=size(y1,1);
        for j=1:n1
           [a,b]=mycompare(x,y1(j,:)); 
           y((j-1)*m+1:j*m,1:i-1)=repmat(y1(j,1:i-1),m,1);
           y((j-1)*m+1:j*m,i)=a;
        end
        y1=y(1:n1*m,1:i);
    end
elseif iscell(x)&&ischar(x{1})    %��ĸ��ʾ
    n=length(x);
    temp=cell(n,1);
    for i=1:n
        temp{i}=1:n;
    end
    y1=exhaus(temp);
    r=size(y1,1);
    num=[];
    for i=1:r
        m=1;
        for j=1:n-1
            for k=j+1:n
                if y1(i,j)==y1(i,k)
                    m=m+1;
                end
            end 
        end
        if m>=2
           num=[num i];
        end
    end
    y=redu(y1,num,'r');
    y1=(1:n)';
    r=size(y,1);
    y2=cell(1,1);
    for i=1:r
        for j=1:n
            a=find(y(i,j)==y1);
            y2(i,j)=x(a);
        end
    end
end   
        
        

