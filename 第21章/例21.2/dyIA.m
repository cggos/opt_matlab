function [best_x,fval]=dyIA(fun,popsize,iter_max,r,LB,UB)   %动态免疫算法,求极小
if isempty(iter_max)
    iter_max=1000;
end
if isempty(popsize)||popsize<100
    popsize=100;
end
if isempty(r)
    r=0.9;
end
[numvar,c]=size(LB);     %多类
for i=1:popsize      %初始化
    if c==1
        pop(i).x=LB'+rand(1,numvar).*(UB-LB)';
    else
        pop(i).x=LB'+rand(c,numvar).*(UB-LB)';
    end
    pop(i).fitness=fun(pop(i).x);
    y1(i)=pop(i).fitness;
end
[a,b1]=sort(y1,'ascend');
cbest=pop(b1(1));   %求极小
sigma=unifrnd(0.1,0.5);
Ag=cbest;
for iter=1:iter_max
    d1=zeros(popsize,popsize);  %距离差
    f=zeros(popsize,popsize);  %适应度差
    for i=1:popsize
        d(i,1)=norm(pop(i).x-Ag.x);
        for j=1:popsize
            if j~=i
               d1(i,j)=norm(pop(i).x-pop(j).x);
               f(i,j)=pop(i).fitness/pop(j).fitness;   %用比值
            end
        end
    end
    for i=1:popsize   %浓度
       temp1=d1(i,:);temp2=f(i,:);
       temp1(i)=inf;temp2(i)=inf;
       a2=find(temp2(find(temp1<r))<1.02);
       n1=length(a2);
       pop(i).c=n1/popsize;
       pop(i).aff=1/(1+d(i,1))*pop(i).fitness/pop(i).c;   %亲和度
       y2(i)=pop(i).aff;
    end    
    [a,b2]=sort(y2,'descend');
    pop=pop(b2);
    clear b2 y2
    I=ceil(unifrnd(0.5,0.6)*popsize);
    for i=1:I                %突变
        alpha=1-exp(-norm(pop(i).x-Ag.x));
        pop(i).x=pop(i).x+unifrnd(0,alpha)*(Ag.x-pop(i).x);
        pop(i).x=boundtest(pop(i).x,LB,UB);
    end
    newpop1=delsample(pop(1:I));
    I1=ceil(unifrnd(0.8,0.9)*I);
    k=randperm(I);
    newpop=newpop1(k(1:I1));   %随机选
    retainpop=pop(I+1:end);
    newpop2=[];
    while 1
        temp=retainpop(1);
        L=size(retainpop,2);
        for i=2:L
            f1(i)=abs(retainpop(i).fitness-temp.fitness); 
        end
        f1(1)=inf;
        m=find(f1(1:L)<sigma);
        if length(m)~=1
            newpop2=[newpop2 temp retainpop(m(1:end-1))];
            retainpop=redu(retainpop,[1 m],'c');
        elseif length(m)==1
            newpop2=[newpop2 temp];
            retainpop=redu(retainpop,[1 m],'c');
        elseif isempty(m)
            newpop2=[newpop2 temp];
            retainpop=redu(retainpop,1,'c');
        end
        if size(retainpop,2)==1||isempty(retainpop)
            break
        end
    end
    pop=[newpop newpop2];
    Gn=size(pop,2);
    num=ceil(popsize*unifrnd(0.05,0.08));
    for i=1:num
        if c==1
           pop(Gn+i).x=LB'+rand(1,numvar).*(UB-LB)';
        else
           pop(Gn+i).x=LB'+rand(c,numvar).*(UB-LB)';
        end
    end
    popsize=Gn+num;
    for i=1:popsize
        pop(i).fitness=fun(pop(i).x);
        y3(i)=pop(i).fitness;
    end
    [a3,b3]=sort(y3,'ascend');
    if a3(1)<cbest.fitness
        cbest=pop(b3(1));
    end
    pop(b3(end))=pop(b3(1));
    temp3=pop;
    clear a3 b3 y3 d1 f k retainpop newpop1 newpop f1 d pop m
    pop=temp3;
    popsize=size(pop,2);
    clear temp3 temp1 temp2 L
    if cbest.fitness<Ag.fitness
       Ag=cbest;
    end
end
best_x=Ag.x;
fval=Ag.fitness;



    

