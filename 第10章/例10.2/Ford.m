function [f,nf,wf,nwf,No]=Ford(C,Flow)   %图论中的Ford算法求解最大流
   %f         ：最大流量矩阵
   %nf        ：最大净流量矩阵
   %wf         ：最大流量
   %nwf        ：最大净流量
   %No                ：标号函数，可得最小割
   %C                 ：弧上的容量
   %Flow              ：弧上的现有流量矩阵，默认情况下为0

n=size(C,2);
if nargin == 1   
    f=zeros(n,n);
else
    f=Flow;
end

%（2）进行标号过程运算
No=zeros(1,n);
d=zeros(1,n);
while 1
    No(1)=n + 1;
    d(1)=inf;         %给源标号
    while 1
        pd=1;         %标号过程
        for i=1:n
            if No(i)    %选择一个已标号的点
                for j=1:n
                    if No(j)==0 && f(i,j)<C(i,j)  %对于未标号的点j，当弧（i,j）非饱和时
                        No(j)=i;
                        d(j)=C(i,j)-f(i,j);
                        pd=0;
                        if d(j)>d(i)
                            d(j)=d(i);
                        end
                    elseif No(j)==0 && f(j,i)>0    %对于未标号的点j，当弧（i,j）非零流时
                        No(j)=-i;
                        d(j)=f(j,i);
                        pd=0;
                        if d(j)>d(i)
                            d(j)=d(i);
                        end
                    end
                end
            end
        end
        if (No(n)||pd)
            break;
        end
    end                 %当汇得到标号或无法标号时，则终止标号过程
    if (pd)
        break;
    end                 %当汇未得到标号，Flow已经为最大流，则算法终止
    %进入流量调整过程，dvt表示调整量
    dvt=d(n);
    t=n;
    while 1
        if No(t)> 0
            f(No(t),t)=f(No(t),t) + dvt;     %前向弧调整
        elseif No(t)<0
            f(-No(t),t)=f(-No(t),t) - dvt;   %后向弧调整
        end
        if No(t)==1
            for i=1:n
                No(i)=0;
                d(i)=0;
            end
            break
        end                   %当t的标号为源时，则终止调整过程
        t=No(t);
    end
end                           %继续调整前一段弧上的流Flow
%（3）求最大净流量矩阵
if nargin==1
   nf=f;
else
    nf=f-Flow;  %最大净流量矩阵
end
%（4）计算最大净流量
wf=0;
nwf=0;
for j=1:n
    wf=wf+f(1,j);
    nwf=nwf+nf(1,j);
end