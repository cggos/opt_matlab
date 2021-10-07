function out=myfzeros(phi,x0)   %phi应为符号表达式
d1=diff(phi);
d2=diff(d1);
ds=char(d1);
if findstr(ds,'*')
   ds=strrep(ds,'*','.*');
end
if findstr(ds,'^')
   ds=strrep(ds,'^','.^');
end
if findstr(ds,'/')
   ds=strrep(ds,'/','./');
end
xs=linspace(x0(1),x0(2),300);
for i=1:300
  fx(i)=feval(inline(ds),xs(i));
end
FZ=fx(1:end-1).*fx(2:end);
Ik=find(FZ<=0);
for k=1:length(Ik)
    y(k)=fzero(inline(ds),xs(Ik(k)));
end
d2f=subs(d2,y);
xM=y(d2f<0);
xm=y(d2f>0);
if exist('xM','var')&&exist('xm','var')
   out.max=xM;
   out.min=xm;
elseif exist('xM','var')&&exist('xm','var')==0
    out.max=xM;
elseif exist('xm','var')&&exist('xM','var')==0
    out.min=xm;
end
out.value=eval(subs(phi,y));
    
