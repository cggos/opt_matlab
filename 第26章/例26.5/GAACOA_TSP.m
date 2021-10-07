function [Shortest_Route,Shortest_Length]=GAACOA_TSP(city,popsize,px,pm,iter_max,Alpha,Beta,Rho,Q)  %遗传－蚁群算法求TSP
[NC,d]=city2d(city);
for i=1:popsize
    pop(i,:)=randperm(NC);
    fitness(i,1)=value(pop(i,:),d);
end
[fval,idex]=min(fitness);
bestx=pop(idex,:);
flag=0;
for iter=1:iter_max(1)
    SumFitValue=sum(fitness);     %所有个体适应值之和
    Ave_x=fitness./SumFitValue;    %所有个体适应值的平均值
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
        PosCutPoint=sort(ceil(NC*rand(1,2)));     %随机确定交叉的切点位置
        while PosCutPoint(1)==PosCutPoint(2)||PosCutPoint(2)-PosCutPoint(2)==NC-1
            PosCutPoint=sort(ceil(NC*rand(1,2)));
        end     
        r1=rand();
        if r1<=px     %进行交叉操作
            nx=pop(FatherSelection,PosCutPoint(1):PosCutPoint(2));
            if PosCutPoint(1)==1
               y2=findpos3(pop(MotherSelection,:),nx);            
               temp=redu(pop(MotherSelection,:),y2,'c');
               new(k,:)=[nx temp];
            else
               new(k,:)=TSPop(pop(FatherSelection,:),pop(MotherSelection,:),'cox');
            end    
            r2=rand();
            if r2<=pm %进行变异操作
                PosMutPoint=ceil(rand*NC);   %随机确定变异元素的位置
                while PosMutPoint==NC
                    PosMutPoint=ceil(rand*NC);
                end
                temp=new(k,PosMutPoint);
                new(k,PosMutPoint)=new(k,PosMutPoint+1);
                new(k,PosMutPoint+1)=temp;
            end
        else
            new(k,:)=pop(FatherSelection,:);
        end
    end
    pop=new;
    for m=1:popsize
        fitness1(m)=value(pop(m,:),d);    %子代个体适应值，
        if iter_max(1)>10
           if iter>1
              r(m)=(fitness1(m)-fitness(m))/fitness(m);
           end
           fitness(m)=fitness1(m);
        end
    end
    [a,b]=min(fitness1);
    if a<fval
        bestx=pop(b,:);
        fval=a;
    end
    if iter_max(1)>10
        if exist('r','var')==1
           R1=r./sum(r);
           if R1<0.05
              flag=flag+1;
           else
              flag=0;
           end
       end
       if flag==10
          break
       end
    end
end
Eta=1./d;          %Eta为启发因子，这里设为距离的倒数
Tau=ones(NC,NC);    %信息素矩阵
Tabu1=bestx;
Delta_Tau= zeros(NC,NC);
for j = 1:(NC-1)
    Delta_Tau(Tabu1(j),Tabu1(j+1)) = Delta_Tau(Tabu1(j),Tabu1(j+1))+Q/fval;
end
Delta_Tau(Tabu1(NC),Tabu1(1)) = Delta_Tau(Tabu1(NC),Tabu1(1))+Q/fval;
Tau=Tau+Delta_Tau;
for i=1:NC-1
    Tau(bestx(i),bestx(i+1))=Tau(bestx(i),bestx(i+1))+Q/fval;
    Tau(bestx(i+1),bestx(i))=Tau(bestx(i),bestx(i+1));
end
Tabu=zeros(popsize,NC);   %存储并记录路径的生成
R_best = zeros(iter_max(2),NC);      %各代最佳路线
L_best = inf.*ones(iter_max(2),1);  %各代最佳路线的长度
L_ave = zeros(iter_max(2),1);       %各代路线的平均长度
for iter=1:iter_max(2)         %停止条件之一：达到最大迭代次数
    Randpos = [];
    for i = 1:(ceil(popsize/NC))
        Randpos = [Randpos,randperm(NC)];
    end
    Tabu(:,1) = (Randpos(1,1:popsize))';
    for j = 2:NC
        for i = 1:popsize
            visited = Tabu(i,1:(j-1));    %已访问的城市
            J = zeros(1,(NC-j+1));         %待访问的城市
            P = J;                        %待访问城市的选择概率分布
            Jc = 1;
            for k = 1:NC
                if isempty(find(visited==k, 1))
                    J(Jc) = k;
                    Jc = Jc+1;
                end
            end
            %下面计算待选城市的概率分布
            for k = 1:length(J)
                P(k) = (Tau(visited(end),J(k))^Alpha)*(Eta(visited(end),J(k))^Beta);
            end
            P=P/(sum(P));
            %按概率原则选取下一个城市
            Pcum = cumsum(P);
            Select = find(Pcum>=rand);
            if isempty(Select)
                to_visit = J(ceil(rand*length(J)));
            else
                to_visit = J(Select(1));
            end
            Tabu(i,j) = to_visit;
        end
    end
    if iter >= 2
        Tabu(1,:) = R_best(iter-1,:); 
    end
    L = zeros(popsize,1);
    for i = 1:popsize
        R = Tabu(i,:);
        for j = 1:(NC-1)
            L(i) = L(i)+d(R(j),R(j+1));
        end
        L(i) = L(i)+d(R(end),R(1));
    end
    [a,b]=sort(L,'ascend');
    temp1=Tabu(b(1),:);
    temp2=Tabu(b(2),:);
    temp3=TSPop(temp1,temp2,'cpmx');
    temp3_F=value(temp3,d);
    if temp3_F<a(1)
        Tabu(b(1),:)=temp3;
        L(b(1),:)=temp3_F;
    end
    for i=1:popsize
        if rand<0.03
           new1=TSPop(Tabu(i,:),'em');
           F1=value(new1,d);
           if F1<L(i)
             L(i)=F1;
             Tabu(i,:)=new1;
           end
        end
    end
    L_best(iter) = min(L);
    pos = find(L==L_best(iter));
    R_best(iter,:) = Tabu(pos(1),:);
    L_ave(iter) = mean(L);
    Delta_Tau = zeros(NC,NC);
    for i = 1:popsize
        for j = 1:(NC-1)
            Delta_Tau(Tabu(i,j),Tabu(i,j+1)) = Delta_Tau(Tabu(i,j),Tabu(i,j+1))+Q/L(i);
        end
        Delta_Tau(Tabu(i,NC),Tabu(i,1)) = Delta_Tau(Tabu(i,NC),Tabu(i,1))+Q/L(i);
    end
    Tau = (1-Rho).*Tau+Delta_Tau;
    Tabu = zeros(popsize,NC);
end
Pos = find(L_best==min(L_best));
Shortest_Route= R_best(Pos(1),:);
Shortest_Length = L_best(Pos(1));
TSPplot(city,Shortest_Route);
