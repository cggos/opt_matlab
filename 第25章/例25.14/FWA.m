function [bestx,fval]=FWA(fun,Nf,Ng,iter_max,LB,UB)     %烟花算法
NC=size(LB,1);
for i=1:Nf    %烟花数
    fireworks(i,:)=LB'+(UB-LB)'.*rand(1,NC);
    f1=fun(fireworks(i,:));
    if f1>=0     %求极小
        fitness(i,1)=1/(f1+1);
    else
        fitness(i,1)=abs(f1)+1;
    end
end
[a1,b1]=sort(fitness,'descend');
best=a1(1);worst=a1(end);
bestx=fireworks(b1(1),:);
fval=a1(1);
M=50;A=40;a=0.04;b=0.8;Aint=0.001;Afinal=0.02;   %算法的参数
for iter=1:iter_max
    total_min=0;total_max=0;
    Amin=Aint-(Aint-Afinal)*sqrt(2*iter_max-iter)/iter_max;
    for i=1:Nf
        total_min=total_min+fitness(i,1)-best;
        total_max=total_max+worst-fitness(i,1);
    end
    for i=1:Nf
        Ai=A*(fitness(i,1)-best+eps)/(total_min+eps);   %半径
        if Ai<Amin
            Ai=Amin;
        end
        s=M*(worst-fitness(i,1)+eps)/(total_max+eps);   %火花数
        if s<a*M
            s=round(a*M);
        elseif s>b*M
            s=round(b*M);
        else
            s=round(s);
        end
        pop=[];y1=[];
        for j=1:s         %烟花爆炸
            new=fireworks(i,:);
            num=ceil(NC*rand);
            z=randperm(NC);
            z=z(1:num);
            for k=1:num
               new(z(k))=new(z(k))+unifrnd(-1,1)*Ai;   %火花的位置
               if new(z(k))>UB(z(k),1)
                  new(z(k))=LB(z(k),1)+rand*(UB(z(k),1)-LB(z(k),1));
               end
               if new(z(k))<LB(z(k),1)
                  new(z(k))=LB(z(k),1)+rand*(UB(z(k),1)-LB(z(k),1));
               end      
            end
            pop=[pop;new];
            f2=fun(new);
            if f2>=0
                f2=1/(f2+1);
            else
                f2=abs(f2)+1;
            end
            y1=[y1;f2];
        end
    end
    for i=1:Ng      %高斯烟花
        n1=ceil(rand*Nf);
        new=fireworks(n1,:);
        num=ceil(NC*rand);
        z=randperm(NC);
        z=z(1:num);
        for k=1:num
            new(z(k))=new(z(k))*randn;
            if new(z(k))>UB(z(k),1)
                new(z(k))=LB(z(k),1)+rand*(UB(z(k),1)-LB(z(k),1));
            end
            if new(z(k))<LB(z(k),1)
                new(z(k))=LB(z(k),1)+rand*(UB(z(k),1)-LB(z(k),1));
            end
        end
        pop=[pop;new];   
        f2=fun(new);
        if f2>=0
           f2=1/(f2+1);
        else
           f2=abs(f2)+1;
        end
        y1=[y1;f2];
    end
    temp=[fireworks;pop];
    y=[fitness;y1];
    [a2,b2]=sort(y,'descend');
    temp=temp(b2,:);
    y=y(b2,1);
    N=size(temp,1);
    m=ceil(0.2*N);
    G=temp(1:m,:);       %精英
    W=temp(N-m+1:end,:);   %差的个体
    Wy=y(N-m+1:end,1);
    fmin=a2(end);
    for i=1:m
        temp1=W(i,:);
        if Wy(i,1)==fmin
            alpha=1;
        else
            alpha=1/(1+2^(Wy(i,1)/fmin))+0.5;
        end
        n2=ceil(rand*m);
        for k=1:NC
            if rand<alpha
                temp1(k)=0.2*temp1(k)+0.8*G(n2,k);  %变异
            end
        end
        f2=fun(temp1);
        if f2>=0
           f2=1/(f2+1);
        else
           f2=abs(f2)+1;
        end
        if f2>Wy(i,1)
           W(i,:)=temp1;
           Wy(i,:)=f2;
        end
    end
    temp=[temp(1:N-m,:);W];
    y=[y(1:N-m,1);Wy];
    [a3,b3]=max(y);
    d=squareform(pdist(temp));  
    totalR=0;
    for i=2:N
        R(i-1)=sum(d(i,:));
        totalR=totalR+R(i-1);
    end
    R=R./totalR;    %概率
    fireworks(1,:)=temp(b3,:);
    fitness(1,1)=a3;
    for i=2:Nf
        p=rand;
        for j=2:N         %轮盘赌策略            
           if p<=sum(R(1:j-1))            
              fireworks(i,:)=temp(j,:);
              fitness(i,1)=y(j,1);
              break
           end        
        end
    end
    [a1,b1]=sort(fitness,'descend');
    best=a1(1);worst=a1(end);
    if a1(1)>fval
       bestx=fireworks(b1(1),:);
       fval=a1(1);
    end
    clear temp y a2 b2 R d z G W Wy 
end
fval=fun(bestx);
    
        
    
        
    
    
        
        