function [xmin,minf,stats]=newls1(fun,x1,b0,x0,x_syms,b_syms,type,esp)   %��˹��ţ����С���˷�
%x_syms,b_syms�ֱ�Ϊ�������ع鳣�������ƣ�����Ϊƽ���͵���ʽ
%x,y�ֱ�Ϊ�Ա������������ʵ������
if nargin==7
    esp=1e-6;
end
n=size(x1,1);
num=size(b_syms,2);
for i=1:n
   f(i,:)=subs(fun,x_syms,[x1(i,:) b0(i)]);    %�����滻
end
if type==1
  jf=myjacobian1(f,b_syms);    %Ŀ�꺯�����ݶ�
  k=0;
  while k<500
        fk=eval(subs(f,b_syms,x0));
        jfk=eval(subs(jf,b_syms,x0));
        dk=-(jfk'*jfk)^-1*jfk'*fk;
        gk=jfk'*fk;
        if norm(gk)<esp
            xmin=x0;
            break;
        end
        beta=0.5;sigma=0.4;
        m=0;mk=0;
        while m<20
            fnew=0.5*norm(eval(subs(f,b_syms,x0+beta^m*dk')))^2;
            fold=0.5*norm(eval(subs(f,b_syms,x0)))^2;
           if fnew<fold+sigma*beta^mk*gk'*dk
              mk=m;
              break;
           end
           m=m+1;
        end
        x0=x0+beta^mk*dk';
        k=k+1;
  end
end
if type==2   %����ֱ�ӽⷽ�̷�
  f1=subs(fun,'b',0);   %���ݺ������ʽ�ı�
  JF=myjacobian1(f1,b_syms);
  kk=0;
  while kk<300
     H=zeros(num,num);
     for i=1:num
        for j=1:num
          if j==i
             ff=subs(JF(1,i),b_syms,x0);
             for k=1:n
                H(i,j)=H(i,j)+subs(ff^2,x_syms(1:end-1),x1(k,:));
             end
          else
             ff=subs(JF(1,i)*JF(1,j),b_syms,x0);
             for k=1:n
                H(i,j)=H(i,j)+subs(ff,x_syms(1:end-1),x1(k,:));
             end
          end
        end
     end
     ff=subs(f1,b_syms,x0);
     R=eval(subs(ff,x_syms(1:end-1),x1))-b0;
     B=zeros(num,1);
     for i=1:num
        ff=subs(JF(1,i),b_syms,x0);
        for j=1:n
           B(i,1)=B(i,1)-subs(ff,x_syms(1:end-1),x1(j,:))*R(j);
        end
     end
     dk=H\B;
     xk=x0+dk';
     if norm(xk-x0)<esp
        xmin=x0;
        break
     end
     x0=xk;
     kk=kk+1;
  end
end
muk=norm(eval(subs(f,b_syms,x0)));
minf=0.5*muk^2;
if nargout==3
  alpha=[0.05 0.01];
  fT=n-1;fA=1;fe=n-2;
  ff=subs(fun,b_syms,xmin);
  for i=1:size(x1,1)
      yhat(i,1)=(eval(subs(ff,x_syms,[x1(i) 0])));
  end
  SSE=(yhat-b0)'*(yhat-b0);
  SST=(b0-mean(b0))'*(b0-mean(b0));
  Sy=sqrt(SSE/fe);
  if SSE>SST
    R2=0;SSR=0;
  else
    R2=1-SSE/SST;
    SSR=SST-SSE;
  end
  Fb=SSR/SSE/fA*fe;
  F=finv(1-alpha,fA,fe);
  if Fb>max(F)
    tst='�߶�����';
  elseif Fb>min(F)&&Fb<=max(F)
    tst='����';
  else
    tst='������';
  end
  table=cell(4,7);
  table(1,:)={'������Դ','ƫ��ƽ����','���ɶ�','����','Fֵ','F��','������'};
  table(2,:)={'�ع�',SSR,fA,SSR/fA,Fb,min(F),tst};
  table(3,:)={'ʣ��',SSE,fe,SSE/fe,[],max(F),[]};
  table(4,1:3)={'�ܺ�',SST,fT};
  stats={table,sqrt(R2),Sy};
end






