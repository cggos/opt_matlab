function y=cluster_dis(x,pattern,m_center)%各个样本到各类中心的距离及类间距离
n=size(x,1);
[centernum,c]=size(m_center);%
d=zeros(n,centernum);
for j=1:centernum
    for i=1:n
        for k=1:c
          d(i,j)=d(i,j)+(x(i,k)-m_center(j,k))^2;  %各个样本到各类中心的距离
        end
    end
end

a1=zeros(1,centernum);  %各类内距离，要小
for j=1:centernum
    for i=1:n
         if pattern(i)==j
           a1(1,j)=a1(1,j)+d(i,pattern(i));
         else
             continue
         end
    end
end
for i=1:centernum
    n1(i)=length(find(pattern==i));
end
for i=1:centernum
    if isempty(n1)
        continue
    else
      a1(1,i)=a1(1,i)/n1(i);
    end
end
d1=sum(a1);

%类间距离,要大
a2=zeros(1,centernum);
for i=1:centernum
    for j=1:centernum
        if i==j
            continue
        else
          for k=1:c
            a2(i)=a2(i)+(m_center(i,k)-m_center(j,k))^2;
          end
        end
    end
end
d2=sum(a2)/2;

%y=d2/d1;
%y=1/sum(a1); 
%y=0.2*max(a1)-0.8*min(a2);
%y=(0.2*max(a1))/(0.6*min(a2));
y=0.3*min(a1)-0.7*max(a2);
%y=0.2*sum(a1)/(0.8*sum(a2));
%y=1/max(a1);


