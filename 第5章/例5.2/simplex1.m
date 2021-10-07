function minx=simplex1(A,b,c)   %修正单纯形法
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
    c=[c 1000*ones(1,n)];
    m=findeye(A);
    [a,b6]=sort(m(:,1));
    m=m(b6,:);
    m=m(end-r1+1:end,:);
end
c1=size(A,2);
minx1=zeros(1,c1);
cb=c(m(:,1));
B=A(:,m(:,1));
E=B;
xb=b(m(:,2));
B1=1\B;
k=0;
while k<200
    p=cb*B1;
    sigma=c-p*A;
    s=find(sigma<-1e-4);
    if isempty(s)
        minx1(m(:,1))=xb;
        minx2=[];
        if ~isempty(n1)
            minx2=redu(minx1,n1,'c');
            c=redu(c,n1,'c');
        end
        if length(c)>length(c_temp)
            n2=length(c_temp)+1:length(c_temp)+n;
            minx2=redu(minx1,n2,'c');
            c=redu(c,n2,'c');  
        end
        if isempty(minx2)
            minf=minx1*c';
        else
            minf=minx2*c';
        end
        break;
    end
    y=B1*A(:,s(1));
    a1=find(y>0);
    if isempty(a1)
        error('无最优解');
    end
    [a,idex]=min(xb(a1)./y(a1));
    zuyan=y(idex);
    E(:,idex)=-y./zuyan;
    E(idex,idex)=1/zuyan;
    B1=E*B1;
    xb=E*xb;
    m(idex,1)=s(1);
    cb(idex)=c(s(1));
    E=B;
    k=k+1;
end
if k==200
    error('陷入死循环');
end
minx=struct;
minx.x=minx1;
minx.x1=minx2;
minx.f=minf;
    


