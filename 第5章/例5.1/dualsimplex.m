function [minx,b]=dualsimplex(A,b,c)     %对偶单纯形法
m=findeye(A);
r1=size(A,1);
c_temp=c;
n1=find(c_temp==0);
if isempty(m)||size(m,1)<r1
    n=r1-size(m,1);
    y5=mycompare(1:r1,m(:,2));   %找单位矩阵
    c2=zeros(r1,n);
    for i=1:n
       c2(y5(i),i)=1;
    end
    A=[A c2];
    c=[c zeros(1,n)];
    m=findeye(A);
    [a,b6]=sort(m(:,1));
    m=m(b6,:);
    m=m(end-r1+1:end,:);
end
[r,c1]=size(A);
minx1=zeros(1,c1);
cb=c(m(:,1));
f0=cb*b;
A_temp=A;
b_temp=b;
f1=f0;
for i=1:c1
   sigma(i)=cb*A_temp(:,i)-c(i);
end
sigma1=sigma;
k1=0;
while k1<200
    s1=find(b_temp<0);
    if isempty(s1)
        break;
    end
    a9=min(b_temp(s1));
    [a9,k]=mycompare1(b_temp(s1),a9);
    k=s1(k);
    b7=find(A_temp(k,:)<0);
    if isempty(b7)
        error('无最优解');
    end
    a6=min(sigma(b7)./A_temp(k,b7));
    [y5,y6]=mycompare1(sigma(b7)./A_temp(k,b7),a6);
    l=b7(min(y6));
    zuyan=A(k,l);
    for i=1:r
       for j=1:c1 
           if i==k
               A(i,j)=A_temp(k,j)/zuyan;
           else           
               A(i,j)=A_temp(i,j)-A_temp(k,j)*A_temp(i,l)/zuyan;
           end
        end
    end
    for i=1:r
        if i==k
            b(i)=b_temp(k)/zuyan;
        else           
            b(i)=b_temp(i)-b_temp(k)*A_temp(i,l)/zuyan;
        end 
    end
    for i=1:c1
        sigma(i)=sigma1(i)-A_temp(k,i)*sigma1(l)/zuyan;
    end
    f0=f1-b_temp(k)*sigma1(l)/zuyan;
    A_temp=A;
    b_temp=b;
    sigma1=sigma;
    f1=f0;
    k1=k1+1;
end
if k1==200
    error('陷入死循环');
end
y=findeye(A);
minx1(y(:,1)')=b(y(:,2));
minx2=[];
if ~isempty(n1)
    minx2=redu(minx1,n1,'c');
end
if length(c)>length(c_temp)
     n2=length(c_temp)+1:length(c_temp)+n;
     minx2=redu(minx1,n2,'c');
end
minf=f0;
minx=struct;
minx.x=minx1;
minx.x1=minx2;
minx.f=minf;
minx.A=[A b;sigma f0];