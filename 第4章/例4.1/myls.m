function [beta,res,stats]=myls(A,b,type)   %������С����
if nargin==2
    type=1;
end
[n,p]=size(A);m=p+1;
n1=size(b,1);
if n1~=n
    b=b(:);
end
answer=questdlg('�г�������','��ѡ��','��','û��','��');
if strcmp(answer,'��')
   A1=[ones(n,1) A];
elseif strcmp(answer,'û��')
    A1=A;
end
if type==1
    [Q,R]=qr(A1);f=Q'*b;
    beta=R\f;
elseif type==2
    [U,S,V]=svd(A);
    r=rank(S);
    beta=zeros(p,1);
    for i=1:r
       beta=beta+(U(:,i)'*b/S(i,i))*V(:,i);
    end
end
res=norm(b-A1*beta);
alpha=[0.05 0.01];
yhat=A1*beta;
SSR=(yhat-mean(b))'*(yhat-mean(b));
SSE=(yhat-b)'*(yhat-b);
SST=(b-mean(b))'*(b-mean(b));
if p==1
   [ab,tm1,r1,rint,stat]=regress(b,A1);
   Fb=SSR/SSE*(n-2);
   Falpha=finv(1-alpha,1,n-2);
   table=cell(4,7);
   table(1,:)={'������Դ','ƫ��ƽ����','���ɶ�','����','Fֵ','F��','������'};
   table(2,1:6)={'�ع�',SSR,1,SSR,Fb,min(Falpha)};
   table(3,1:6)={'ʣ��',SSE,n-2,SSE/(n-2),[],max(Falpha)};
   table(4,1:3)={'�ܺ�',SST,n-1};
   if Fb>=max(Falpha)
      table{2,7}='�߶�����';
   elseif (Fb<max(Falpha))&&(Fb>=min(Falpha))
       table{2,7}='����';
   else
       table{2,7}='������';
   end
   Sy=sqrt(table{3,4});
   r2=SSR/SST;
   stats={table,sqrt(r2),Sy,};
elseif p>=2
    Fb=SSR/(m-1)/SSE*(n-m);
    Falpha=finv(1-alpha,m-1,n-m);
    table=cell(p+4,7);
    table(1,:)={'������Դ','ƫ��ƽ����','���ɶ�','����','Fֵ','F��','������'};
    table(2+p,1:6)={'�ع�',SSR,m-1,SSR/(m-1),Fb,min(Falpha)};
    table(3+p,1:6)={'ʣ��',SSE,n-m,SSE/(n-m),[],max(Falpha)};
    table(4+p,1:3)={'�ܺ�',SST,n-1};
    if Fb>=max(Falpha)
       table{2+p,7}='�߶�����';
    elseif (Fb<max(Falpha))&&(Fb>=min(Falpha))
       table{2+p,7}='����';
    else
       table{2+p,7}='������';
    end
    R2=SSR/SST;
    R=sqrt(R2);
    Sy=sqrt(SSE/(n-m));
    mnX=mean(b);
    MNX=repmat(mnX,n,1);
    Ljj=diag((b-MNX)'*(b-MNX));
    pj=abs(beta(2:end).*sqrt(Ljj/SST));  %��׼ƫ�ع�ϵ��
    C=diag(inv(A1'*A1));bj2=beta.*beta;
    ssj=bj2(2:end)./C(2:end);
    Fj=ssj/SSE*(n-m);
    Falpha=finv(1-[0.05,0.01],1,n-m);
    ind2=find(Fj>=Falpha(2));
    ind1=find((Fj>=Falpha(1))&(Fj<Falpha(2)));
    ind0=find(Fj<Falpha(1));
    xxx=zeros(size(Fj));
    xxx(ind2)=2;
    xxx(ind1)=1;
    [tem,zbx]=min(Fj);
    xzh={'������','����','�߶�����'};
    for kk=1:p
        table(kk+1,:)={['x',num2str(kk)],ssj(kk),1,ssj(kk),Fj(kk),[],xzh{1+xxx(kk)}};
    end
    table{2,6}=Falpha(1);table{3,6}=Falpha(2);
    stats={table,R,Sy,pj};
end





