function [y,fval,rd2]=dyprog(x,d,decisfun,subfun,trafun,objfun,str)
%x               状态变量, 一列代表一个阶段状态，根据具体问题分析后得出X的初始值取值范围；
%decisfun(k, X)  由阶段k的状态变量x求出相应的允许决策变量；
%subfun       阶段指标函数；
%trafun(k,X,u) 状态转移函数, 其中X是阶段k的某状态变量, u是相应的决策变量；
%objfun(v,f)     第k阶段至最后阶段指标函数, 当ObjFun(v,f)=v+f时, 输入ObjFun可以省略；否则需要给出定义。
%y           结果输出由4列构成,y=[序号组；最优策略组；最优轨线组；指标函数值组]；
%fval            一个列向量, 各元素分别表示y各最优策略组对应始端状态X的最优函数值。
if nargin==5
    str=[];
    objfun=[];
end
if nargin==6
    str=[];
end
k=size(x,2);               %计算动态规划的阶段数，即为数组X的列数
x_isnan=~isnan(x);              %给x的各元素赋初值
t_vub=inf; 
t_vubm=inf*ones(size(x));
f_opt=nan*ones(size(x));        %nan表示无意义元素
d_opt=f_opt;
tmp1=find(x_isnan(:,k));
tmp2=length(tmp1);
for i=1:tmp2
    u=decisfun(k,x(i,k)); %执行decisfun
    tmp3=length(u);
    for j=1:tmp3
        tmp=subfun(k,x(tmp1(i),k),u(j),d);  %执行subfun
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
            tmp40=trafun(ii,x(tmp10(i),ii),u(j));  %执行trafun
            tmp50=x(:,ii+1)-tmp40;
            tmp60=find(tmp50==0);
            if ~isempty(tmp60)
                if isempty(objfun)
                    tmp00=tmp00+f_opt(tmp60(1),ii+1);
                else
                    tmp00=objfun(tmp00,f_opt(tmp60(1),ii+1)); %执行objfun
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
    tmpf(i)=subfun(1,tmpX(i),tmpd(i),d);                   %执行SubObjFun
    y(k*(i-1)+1,[1,2,3,4])=[1,tmpX(i),tmpd(i),tmpf(i)];
    for ii=2:k
        tmpX(i)=trafun(ii-1,tmpX(i),tmpd(i));             %执行TransFun
        tmp1=x(:,ii)-tmpX(i);
        tmp2=find(tmp1==0);
        if ~isempty(tmp2)
            tmpd(i)=d_opt(tmp2(1),ii);
        end
        tmpf(i)=subfun(ii,tmpX(i),tmpd(i),d);              %执行SubObjFun
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
              rd2=[rd2 syms(y2) '→'];
           else
              rd2=[rd2 syms{y2} '→'];
           end
       end
    end
end

