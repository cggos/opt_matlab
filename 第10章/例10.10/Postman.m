function [Edge,Sum]=Postman(M)     %邮递员问题
%Edge      ：边的集合
%Sum       ：行走路径的总长度
%M         ：图的权值矩阵
n=length(M);
B=M;
B(B == inf)=0;                   %如果为无穷大，则赋值为0
B(B~=0)=1;                     %如果不为零，则赋值为1
m=0;
a=sum(B);
edges=sum(sum(B))/2;             %计算图的边数
ed=zeros(2,edges);               %给存放Euler回路的矩阵赋初值
vexs=zeros(1,edges+1);
Matr=B;
for i=1:n
    if mod(a(i),2)==1            %查找次为奇数的顶点（奇点）
        m=m+1;
    end
end

if m~=0
    fprintf('图中奇点的个数为%i个：\n',m);
    fprintf('不存在Euler环游路径.\n');
    Edge = 0;
    Sum = 0;
end

if m==0
    vet=1;
    flag=0;
    R=find(Matr(vet,:)==1);
    for p=1:length(R)
        ed(:,1)=[vet,R(p)];
        vexs(1,1)=vet;
        vexs(1,2)=R(p);
        Matr(vexs(1,2),vexs(1,1))=0;
        flagg=1;
        temp=1;
        while flagg
            [flagg,ed]=Fleury_edf(Matr,edges,vexs,ed,temp);   %调用Fleury_edf函数
            temp=temp+1;
            if ed(1,edges)~=0 && ed(2,edges)~=0
                Edge=ed;
                Edge(2,edges)=1;
                Sum=0;
                for g=1:edges
                    Sum=Sum+M(Edge(1,g),Edge(2,g));
                end
                flagg=0;
                break;
            end
        end
    end
end

function [flag,ed]=Fleury_edf(Matr,edges,vexs,ed,temp)
flag = 1;
for i = 2:edges
    [dvex,f]=Fleury_flecvexf(Matr,i,vexs,edges,ed,temp);   %调用Ch12_Fleury_flecvexf函数
    if f==1
        flag=0;
        break;
    end
    if dvex~=0
        ed(:,i)=[vexs(1,i) dvex];
        vexs(1,i+1)=dvex;
        Matr(vexs(1,i+1),vexs(1,i))=0;
    else
        break;
    end
end

function [dvex,f]=Fleury_flecvexf(Matr,i,vexs,edges,ed,temp)
f=0;
edd=find(Matr(vexs(1,i),:) == 1);
dvex=0;
dvex1=[];
ded=[];
if length(edd)==1
    dvex=edd;
else
    dd=1;
    dd1=0;
    kkk=0;
    for kk=1:length(edd)
        m1=find(vexs==edd(kk));
        if sum(m1)==0
            dvex1(dd)=edd(kk);
            dd=dd+1;
            dd1=1;
        else
            kkk=kkk+1;
        end
    end
    if kkk==length(edd)
        tem=vexs(1,i)*ones(1,kkk);
        edd1=[tem;edd];
        for q1=1:kkk
            w=0;
            ddd=1;
            for q2=1:edges
                if edd1(1:2,q1)==ed(1:2,q2);
                    w=w+1;
                end
            end
            if w==0
                ded(ddd)=edd(q1);
                ddd=ddd+1;
            end
        end
    end
    
    if temp<=length(dvex1)
        dvex=dvex1(temp);
    elseif temp>length(dvex1) && temp<=length(ded)
        dvex=ded(temp);
    else
        f=1;
    end
end