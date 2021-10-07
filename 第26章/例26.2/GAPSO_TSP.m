function [MinDistance,Path] = GAPSO_TSP(city,popsize,MaxGen)
[n,d]=city2d(city); %n表示问题的规模（城市个数）     
individual=zeros(popsize,n);
for i=1:popsize
    individual(i,:)=randperm(n);    %随机产生粒子位置
    indiFit(i,1)=value(individual(i,:),d);
end
[value1,index]=min(indiFit);
tourPbest=individual;                                    %当前个体最优
tourGbest=individual(index,:) ;                          %当前全局最优
xnew1=individual;
trace=zeros(1,MaxGen);   %存放各代最优路径的长度及各代路径的平均长度
for N = 1:MaxGen
    for i = 1:popsize
       % 与个体最优进行交叉
        c1 = unidrnd(n-1);    %产生交叉位
        c2 = unidrnd(n-1);    %产生交叉位
        while c1 == c2
            c1 = round(rand*(n-2))+1;
            c2 = round(rand*(n-2))+1;
        end
        chb1 = min(c1,c2);
        chb2 = max(c1,c2);
        cros = tourPbest(i,chb1:chb2);
        ncros = size(cros,2);      
        %删除与交叉区域相同元素
        for j = 1:ncros
            for k = 1:n
                if xnew1(i,k) == cros(j)
                    xnew1(i,k) = 0;
                    for t = 1:n-k
                        temp = xnew1(i,k+t-1);
                        xnew1(i,k+t-1) = xnew1(i,k+t);
                        xnew1(i,k+t) = temp;
                    end
                end
            end
        end
        %插入交叉区域
        xnew1(i,n-ncros+1:n) = cros;
        %新路径长度变短则接受
        Dist = 0;
        for j=1:n-1
            Dist = Dist+d(xnew1(i,j),xnew1(i,j+1));
        end
        Dist = Dist+d(xnew1(i,1),xnew1(i,n));
        if indiFit(i) > Dist
            individual(i,:)=xnew1(i,:);
            indiFit(i,1)=Dist;
        end
        
        % 与全体最优进行交叉
        c1 = round(rand*(n-2))+1;  %产生交叉位
        c2 = round(rand*(n-2))+1;  %产生交叉位
        while c1 == c2
            c1 = round(rand*(n-2))+1;
            c2 = round(rand*(n-2))+1;
        end
        chb1 = min(c1,c2);
        chb2 = max(c1,c2);
        cros = tourGbest(chb1:chb2); 
        ncros = size(cros,2);      
        %删除与交叉区域相同元素
        for j = 1:ncros
            for k = 1:n
                if xnew1(i,k) == cros(j)
                    xnew1(i,k) = 0;
                    for t = 1:n-k
                        temp = xnew1(i,k+t-1);
                        xnew1(i,k+t-1) = xnew1(i,k+t);
                        xnew1(i,k+t) = temp;
                    end
                end
            end
        end
        %插入交叉区域
        xnew1(i,n-ncros+1:n) = cros;
        %新路径长度变短则接受
        Dist = 0;
        for j = 1:n-1
            Dist = Dist+d(xnew1(i,j),xnew1(i,j+1));
        end
        Dist = Dist+d(xnew1(i,1),xnew1(i,n));
        if indiFit(i) > Dist
            individual(i,:) = xnew1(i,:);
            indiFit(i,1)=Dist;
        end
        
       % 变异操作
        c1 = round(rand*(n-1))+1;   %产生变异位
        c2 = round(rand*(n-1))+1;   %产生变异位
        while c1 == c2
            c1 = round(rand*(n-2))+1;
            c2 = round(rand*(n-2))+1;
        end
        temp = xnew1(i,c1);
        xnew1(i,c1) = xnew1(i,c2);
        xnew1(i,c2) = temp;
        
        %新路径长度变短则接受
        Dist = 0;
        for j = 1:n-1
            Dist = Dist + d(xnew1(i,j),xnew1(i,j+1));
        end
        Dist = Dist + d(xnew1(i,1),xnew1(i,n));      %计算环游路径的长度
        if indiFit(i) > Dist
            individual(i,:) = xnew1(i,:);
            indiFit(i,1)=Dist;
        end
    end
    [value2,index] = min(indiFit);
    if value2<value1
       trace(2,N) = indiFit(index);    %最优路径的长度
       tourGbest = individual(index,:); 
       value1=value2;
    end  
end
MinDistance = value1;
Path = tourGbest;
TSPplot(city,Path);