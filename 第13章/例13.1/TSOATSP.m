function cBest=TSOATSP(city,StopL,Candidate_Num,Tabu_Length)
%city为城市坐标，itermax最大迭代次数,L禁忌表长度
[m,Distance_List]=city2d(city);
%m=size(city,1);
PlaceNum=m;                   %旅行地的数量
if nargin==1
    StopL=3000;
    Tabu_Length=60;
    Candidate_Num=40;
elseif nargin==2
    Tabu_Length=60;
    Candidate_Num=40;
elseif nargin==3
    Tabu_Length=60;
end
TabuList=zeros(Tabu_Length,PlaceNum);   	%禁忌表(tabu list)
m_solution.route=randperm(PlaceNum);      	%随机生成初始解序列;
m_solution.fitness=value(m_solution.route,Distance_List);
m_solution.idex=1;
cBest=m_solution;               %最优解
for iterm=1:StopL
    A=zeros(Candidate_Num,2);   
    A=ceil(PlaceNum*rand(Candidate_Num,2));  	    %随机产生元素值取值范围为[0,PlaceNum]的2维行向量，生成候选解
    for i=1:Candidate_Num             %候选解
       for j=1:Candidate_Num
          if j~=i
            while isequal(A(i,:),A(j,:))     
                A(j,:)=ceil(PlaceNum*rand(1,2));  
            end
          end
       end
    end
    S=m_solution.route;
    Si=zeros(Candidate_Num,PlaceNum);
    for i=1:Candidate_Num
        Si(i,:)=S;                                  %设置候选解的每一行S为初始随机序列
        Si(i,[A(i,1),A(i,2)])=S([A(i,2),A(i,1)]);   %交换第i行的两个位置的数据位置的编号由A产生
        %计算第i个候选解的总距离值
        Total_Distance(i)=value(Si(i,:),Distance_List);
    end
    [a,b]=sort(Total_Distance);
    for i=1:Candidate_Num
        Candidate(i).fitness=Total_Distance(b(i));
        Candidate(i).route=Si(b(i),:);
        if iterm==1
            Candidate(i).bTabu=0;
        end 
    end
    for i=1:Candidate_Num
        for j=1:Tabu_Length
            if isequal(Candidate(i).route,TabuList(j,:))
                Candidate(i).bTabu=1;     %在禁忌表中
                break
            else
                Candidate(i).bTabu=0;
            end
        end
    end
    balltabu=1;    %判断是否都在禁忌表中
    for i=1:Candidate_Num
       if Candidate(i).bTabu==0     %不在禁忌表中
           balltabu=0;
           break
       end
    end
    bestCandidate=Candidate(1);   %最优候选解
    if bestCandidate.fitness<cBest.fitness
        m_solution=bestCandidate;
        for i=1:Tabu_Length-1
            TabuList(i,:)=TabuList(i+1,:);
        end
        TabuList(Tabu_Length,:)=bestCandidate.route;
    else
        if balltabu==0
            for k=1:Candidate_Num
                if Candidate(k).bTabu==0
                    m_solution=Candidate(k);
                    for j=1:Tabu_Length-1
                       TabuList(j,:)=TabuList(j+1,:);
                    end
                    TabuList(Tabu_Length,:)=bestCandidate.route;
                    break
               end
            end
        else
            m_solution=bestCandidate;
            for j=1:Tabu_Length-1
                 TabuList(j,:)=TabuList(j+1,:);
            end
            TabuList(Tabu_Length,:)=bestCandidate.route;
        end
    end
    if m_solution.fitness<cBest.fitness
        cBest=m_solution;
        cBest.index=iterm;
    end
end
Route=cBest.route;        %搜索到的最优解的环游路径
TSPplot(city,Route);      %画图
%gname;



