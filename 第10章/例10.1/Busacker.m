function [f,wf,zwf]=Busacker(C,b)    %图论中的Busacker_Gowan迭代算法
%====输出参数====
   %f     ：最小费用最大流矩阵
   %wf              ：最大流量
   %zwf              ：最小费用
%====输入参数====
   %C                    ：弧容量
   %b                    ：弧上单位流量的费用
n=size(C,2);                    %图的节点个数
wf=0;wf0=Inf;      %最大流量MaxFlow及预定的流量值MaxFlow0
f=zeros(n,n);    %取初始可行流为“零流”
while 1
    a=ones(n,n)*inf;
    for i=1:n
        a(i,i)=0;
    end
    %第1步：构造有向赋权图
    for i=1:n
        for j=1:n
            if(C(i,j)>0 &&f(i,j)==0)
                a(i,j)=b(i,j);
            elseif(C(i,j)>0 && f(i,j)==C(i,j))
                a(j,i)=-b(i,j);
            elseif(C(i,j)>0)
                a(i,j)=b(i,j);
                a(j,i)=-b(i,j);
            end
        end
    end
    for i=2:n
        p(i)=inf;
        s(i)=i;
    end
    %第2步：求最短路——求有向赋权图中起点到终点的最短路
    for k=1:n
        pd=1;      
        for i=2:n
            for j=1:n
                if(p(i)>p(j)+a(j,i))
                    p(i)=p(j)+a(j,i);
                    s(i)=j;
                    pd=0;
                end
            end
        end
        if(pd)
            break;
        end;
    end
    
    if(p(n)==inf)
        break;
    end          %不存在起点到终点的最短路, 算法终止。
    dvt=inf;
    t=n;
    %第3步：进入调整过程,dvt表示调整量
    while 1    %计算调整量
        if(a(s(t),t)>=0)
            dvtt = C(s(t),t)-f(s(t),t); %前向弧调整量
        elseif(a(s(t),t)<0)
            dvtt=f(t,s(t));
        end                                            %后向弧调整量
        if(dvt>dvtt)
            dvt=dvtt;
        end
        if(s(t)==1)
            break;
        end      %当t的标号为起点时, 终止计算调整量
        t=s(t);
    end
    %继续调整前一段弧上的流f
    pd=0;
    if(wf+dvt>=wf0)
        dvt=wf0-wf;
        pd=1;
    end
    t=n;
    while 1
        if(a(s(t),t)>=0)
            f(s(t),t)=f(s(t),t)+dvt;     %前向弧调整
        elseif(a(s(t),t) < 0)
            f(t,s(t))=f(t,s(t))-dvt;
        end                                                              %后向弧调整
        if(s(t)==1)
            break;
        end   %当t的标号为起点时, 终止调整过程
        t=s(t);
    end
    if(pd)
        break;
    end       %如果最大流量达到预定的流量值
    %计算最大流量
    wf=0;
    for j=1:n
        wf=wf+f(1,j);
    end;
end
%计算最小费用
zwf=0;
for i=1:n
    for j=1:n
        zwf=zwf+b(i,j)*f(i,j);
    end
end