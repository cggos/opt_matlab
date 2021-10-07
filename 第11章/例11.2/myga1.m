function [max_x,maxfval]=myga1(fun,LB,UB,popsize,iterm_max,px,pm,eps1)
  %fun       ：待优化的目标函数
  %LB            ：自变量下界
  %UB            ：自变量上界
  %popsize       ：种群大小
  %iterm_max      ：最大迭代数
  %px            ：杂交概率
  %pm            ：变异概率
  %eps1           ：自变量离散精度
if isempty(pm)
    pm=0.1;
end
if isempty(px)
    px=0.90;
end
if isempty(iterm_max)
    iterm_max=8000;
end
if isempty(popsize)
    popsize=50;
end
if isempty(LB)&&isempty(UB)
    nvar=input('请输入变量数目nvar＝ ');
else
    nvar=size(LB,1);     %变量数
end
CodeLen=nvar*max(ceil(log2((UB-LB)/eps1 + 1)));      %根据自变量离散精度，确定二进制编码位串的长度
x=zeros(popsize,CodeLen);                      %种群编码的初始值
for i=1:popsize
    x(i,:)=Initial(CodeLen);              %调用子程序Initial初始化种群
    FitValue(i)=fun(Dec1(LB,UB,x(i,:),CodeLen));    %个体适应值
end
for i=1:iterm_max
    SumFitValue=sum(FitValue);     %所有个体适应值之和
    Ave_x=FitValue/SumFitValue;    %所有个体适应值的平均值
    Prob_Ave_x=0;
    Prob_Ave_x(1)=Ave_x(1);
    for j=2:popsize                    %用于轮盘赌策略的概率累加
        Prob_Ave_x(j)=Prob_Ave_x(j-1)+Ave_x(j);
    end
    for k=1:popsize
        sita=rand();
        for n=1:popsize
            if sita<=Prob_Ave_x(n)
                FatherSelection=n;   %根据轮盘赌策略选择父代
                break;
            end
        end
        MotherSelection=floor(rand()*(popsize-1))+1;      %随机确定母亲
        PosCutPoint=floor(rand()*(CodeLen-2))+1;     %随机确定单点交叉的切点位置
        r1=rand();
        if r1<=px     %进行交叉操作
            nx(k,1:PosCutPoint)=x(FatherSelection,1:PosCutPoint);
            nx(k,(PosCutPoint+1):CodeLen)=x(MotherSelection,(PosCutPoint+1):CodeLen);
            r2=rand();
            if r2<=pm %进行变异操作
                PosMutPoint=round(rand()*(CodeLen-1)+1);   %随机确定变异元素的位置
                nx(k,PosMutPoint)=~nx(k,PosMutPoint);
            end
        else
            nx(k,:)=x(FatherSelection,:);
        end
    end
    x=nx;
    for m=1:popsize
        FitValue(m)=fun(Dec1(LB,UB,x(m,:),CodeLen));    %子代个体适应值，调用子函数Dec
    end
    [a,b]=max(FitValue);
    best=x(b,:);
    x(popsize,:)=best;
end
max_x=Dec1(LB,UB,best,CodeLen);
maxfval=fun(max_x);


function result=Initial(length)           %初始化函数
for i=1:length
    r=rand();
    result(i)=round(r);
end


function y=Dec1(LB,UB,x,CodeLen)           %二进制转换为十进制编码
nvar=size(LB,1);
sublen=CodeLen/nvar;
pow_two=2.^(0:sublen)';
maxintval=((2^sublen))-1;
range=UB'-LB';
start=1;
fin=sublen;
for j=1:nvar
    tvars(1:sublen)=x(start:fin);
    start=start+sublen;
    fin=fin+sublen;
    temp1=0;
    for k=1:sublen
         temp1=temp1+pow_two(k)*tvars(sublen-k+1);
    end
    y(j)=(range(j)*(temp1/maxintval))+LB(j);
end