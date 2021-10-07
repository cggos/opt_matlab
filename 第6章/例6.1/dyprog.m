function [y,fval,rd2]=dyprog(x,d,decisfun,subfun,trafun,objfun,str)
%x               ״̬����, һ�д���һ���׶�״̬�����ݾ������������ó�X�ĳ�ʼֵȡֵ��Χ��
%decisfun(k, X)  �ɽ׶�k��״̬����x�����Ӧ��������߱�����
%subfun       �׶�ָ�꺯����
%trafun(k,X,u) ״̬ת�ƺ���, ����X�ǽ׶�k��ĳ״̬����, u����Ӧ�ľ��߱�����
%objfun(v,f)     ��k�׶������׶�ָ�꺯��, ��ObjFun(v,f)=v+fʱ, ����ObjFun����ʡ�ԣ�������Ҫ�������塣
%y           ��������4�й���,y=[����飻���Ų����飻���Ź����飻ָ�꺯��ֵ��]��
%fval            һ��������, ��Ԫ�طֱ��ʾy�����Ų������Ӧʼ��״̬X�����ź���ֵ��
if nargin==5
    str=[];
    objfun=[];
end
if nargin==6
    str=[];
end
k=size(x,2);               %���㶯̬�滮�Ľ׶�������Ϊ����X������
x_isnan=~isnan(x);              %��x�ĸ�Ԫ�ظ���ֵ
t_vub=inf; 
t_vubm=inf*ones(size(x));
f_opt=nan*ones(size(x));        %nan��ʾ������Ԫ��
d_opt=f_opt;
tmp1=find(x_isnan(:,k));
tmp2=length(tmp1);
for i=1:tmp2
    u=decisfun(k,x(i,k)); %ִ��decisfun
    tmp3=length(u);
    for j=1:tmp3
        tmp=subfun(k,x(tmp1(i),k),u(j),d);  %ִ��subfun
        if tmp<=t_vub
            f_opt(i,k)=tmp;
            d_opt(i,k)=u(j);
            t_vub=tmp;
        end
    end
end
for ii=k-1:-1:1
    tmp10=find(x_isnan(:,ii));
    tmp20=length(tmp10);
    for i=1:tmp20
        u=decisfun(ii,x(i,ii));
        tmp30=length(u);
        for j=1:tmp30
            tmp00=subfun(ii,x(tmp10(i),ii),u(j),d); 
            tmp40=trafun(ii,x(tmp10(i),ii),u(j));  %ִ��trafun
            tmp50=x(:,ii+1)-tmp40;
            tmp60=find(tmp50==0);
            if ~isempty(tmp60)
                if isempty(objfun)
                    tmp00=tmp00+f_opt(tmp60(1),ii+1);
                else
                    tmp00=objfun(tmp00,f_opt(tmp60(1),ii+1)); %ִ��objfun
                end
                if tmp00 <= t_vubm(i,ii)
                    f_opt(i,ii)=tmp00;
                    d_opt(i,ii)=u(j);
                    t_vubm(i,ii)=tmp00;
                end
            end
        end
    end
end
fval=f_opt(tmp1,1);
fval=fval(~isnan(fval),1);
y=[];tmpX=[];tmpd=[];tmpf=[];
tmp0=find(x_isnan(:,1));
tmp01=length(tmp0);
for i=1:tmp01
    tmpd(i)=d_opt(tmp0(i),1);
    tmpX(i)=x(tmp0(i),1);
    tmpf(i)=subfun(1,tmpX(i),tmpd(i),d);                   %ִ��SubObjFun
    y(k*(i-1)+1,[1,2,3,4])=[1,tmpX(i),tmpd(i),tmpf(i)];
    for ii=2:k
        tmpX(i)=trafun(ii-1,tmpX(i),tmpd(i));             %ִ��TransFun
        tmp1=x(:,ii)-tmpX(i);
        tmp2=find(tmp1==0);
        if ~isempty(tmp2)
            tmpd(i)=d_opt(tmp2(1),ii);
        end
        tmpf(i)=subfun(ii,tmpX(i),tmpd(i),d);              %ִ��SubObjFun
        y(k*(i-1)+ii,[1,2,3,4]) = [ii,tmpX(i),tmpd(i),tmpf(i)];
    end
end
if ~isempty(str)
    r3=size(str,1);
    rd2=[];
    for j=1:r3
       [syms,n]=findletter(str{j},2);
       [y1,y2]=mycompare1(y(j,2),x(:,j));
        if j==r3
           if n==1
               rd2=[rd2 syms(y2)];
           else
               rd2=[rd2 syms{y2}];
           end
       else
           if n==1
              rd2=[rd2 syms(y2) '��'];
           else
              rd2=[rd2 syms{y2} '��'];
           end
       end
    end
end

