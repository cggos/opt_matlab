function [xmin,minf,mu,lm]=newquad(H,c,Ae,be,Ai,bi,x0,esp)   %光滑牛顿法求解二次规划
n=length(c);
if nargin==6
    esp=1e-6;
    x0=ones(n,1);
end
be=-be;bi=-bi;   %化成标准式
l=length(be);m=length(bi);
beta=0.5;sigma=0.2;gamma=0.05;
ep=0.05;mu=0.05*ones(l,1);
lm=0.05*zeros(m,1);u=[ep;zeros(n+l+m,1)];
k=0;
while k<200
   dh=dah(ep,x0,mu,lm,c,H,Ae,be,Ai,bi);
   if norm(dh)<esp
      xmin=x0;
      break
   end
   A=jacoH(ep,x0,lm,c,H,Ae,be,Ai,bi);
   b=ps(ep,x0,mu,lm,c,H,Ae,be,Ai,bi,gamma)*u-dh;
   dz=A\b;
   de=dz(1);dd=dz(2:n+1);
   if l>0&&m>0
       du=dz(n+2:n+l+1);dl=dz(n+l+2:n+l+m+1);
   end
   if l==0
       dl=dz(n+2:n+m+1);
   end
   if m==0
       du=dz(n+2:n+l+1);
   end
   i=0;
   while i<=20
       t1=beta^i;
       if l>0&&m>0
          dh1=dah(ep+t1*de,x0+t1*dd,mu+t1*du,lm+t1*dl,c,H,Ae,be,Ai,bi);
       end
       if l==0
          dh1=dah(ep+t1*de,x0+t1*dd,mu,lm+t1*dl,c,H,Ae,be,Ai,bi);
       end
       if m==0
           if isempty(dd)
               bb=x0;
           else
               bb=x0+t1*dd;
           end
           dh1=dah(ep+t1*de,bb,mu+t1*du,lm,c,H,Ae,be,Ai,bi);
       end
       if norm(dh1)<=(1-sigma*(1-gamma*ep)*beta^i)*norm(dh)
           mk=i;
           break
       end
       i=i+1;
       if i==20
           mk=10;
       end
   end
   alpha=beta^mk;
   ep=ep+alpha*de;xk=x0+alpha*dd;
   if l>0&&m>0
      mu=mu+alpha*du;lm=lm+alpha*dl;
   end
   if l==0
       lm=lm+alpha*dl;
   end
   if m==0
      mu=mu+alpha*du;
   end
   x0=xk;
   k=k+1;
end
minf=0.5*xmin'*H*xmin+c'*xmin;


function dh=dah(ep,x0,mu,lm,c,H,Ae,be,Ai,bi)
n=length(c);l=length(be);m=length(bi);
dh=zeros(n+l+m+1,1);
dh(1)=ep;
if l>0&&m>0
    dh(2:n+1)=H*x0-Ae'*mu-Ai'*lm+c;
    dh(n+2:n+l+1)=be+Ae*x0;
    for i=1:m
        dh(n+l+1+i)=ph(ep,lm(i),bi(i)+Ai(i,:)*x0);
    end
end
if l==0
    dh(2:n+1)=H*x0-Ai'*lm+c;
    for i=1:m
        dh(n+1+i)=ph(ep,lm(i),bi(i)+Ai(i,:)*x0);
    end
end
if m==0
    dh(2:n+1)=H*x0-Ae'*mu+c;
    dh(n+2:n+l+1)=be+Ae*x0;
end
dh=dh(:);

function y=ph(ep,a,b)
y=a+b-sqrt(a^2+b^2+2*ep^2);


function y=ps(ep,x0,mu,lm,c,H,Ae,be,Ai,bi,gamma)
dh=dah(ep,x0,mu,lm,c,H,Ae,be,Ai,bi);
y=gamma*norm(dh)*min(1,norm(dh));


function A=jacoH(ep,x0,lm,c,H,Ae,be,Ai,bi)
n=length(c);l=length(be);m=length(bi);
A=zeros(n+l+m+1,n+l+m+1);
[dd1,dd2,v1]=ddv(ep,x0,lm,Ai,bi);
if l>0&&m>0
   A=[1,       zeros(1,n),  zeros(1,l), zeros(1,m);
     zeros(n,1),  H,          -Ae',      -Ai';
     zeros(l,1),  Ae,       zeros(l,l), zeros(l,m);
       v1,      dd2*Ai,     zeros(m,l), dd1];
end
if l==0
   A=[1,       zeros(1,n),  zeros(1,m);
     zeros(n,1),  H,         -Ai';
      v1,       dd2*Ai,       dd1];
end
if m==0
   A=[1,        zeros(1,n),  zeros(1,l);
      zeros(n,1),  H,         -Ae';
      zeros(l,1),  Ae,      zeros(l,l)];
end

function [dd1,dd2,v1]=ddv(ep,x0,lm,Ai,bi)
m=length(bi);
dd1=zeros(m,m);dd2=zeros(m,m);v1=zeros(m,1);
for i=1:m
    fm=sqrt(lm(i)^2+(bi(i)+Ai(i,:)*x0)^2+2*ep^2);
    dd1(i,i)=1-lm(i)/fm;
    dd2(i,i)=1-(bi(i)+Ai(i,:)*x0)/fm;
    v1(i)=-2*ep/fm;
end





   
       
       



