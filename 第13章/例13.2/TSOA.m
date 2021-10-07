function [cBest,iterm]=TSOA(fun,varnum,LB,UB,iterm_max,Candidate_Num,L,pm)   %禁忌算法求解函数优化
%fun为函数，itermax最大迭代次数,L禁忌表长度，
%varnum变量数，LB,UB分别为各变量的上下界,pm为变异概率
if nargin==4
    StopL=2000;
    Tabu_Length=10;
    Candidate_Num=300;
    pm=0.75;
elseif nargin==5
    StopL=iterm_max;         	%最大迭代次数
    Tabu_Length=10;
    Candidate_Num=40;
    pm=0.75;
elseif nargin==6
    StopL=iterm_max;
    Tabu_Length=10;
    pm=0.75;
elseif nargin==7
    StopL=iterm_max;
    Tabu_Length=L;
    pm=0.75;
else
    StopL=iterm_max;
    Tabu_Length=L;
end
TabuList=zeros(Tabu_Length,varnum+1);   	%禁忌表(tabu list)
m_solution.idex=1;
m_solution.x=LB'+rand(1,varnum).*(UB-LB)';
m_solution.fitness=fun(m_solution.x);
cBest=m_solution;
for iterm=1:StopL
    if iterm<0.95*StopL
        x0=m_solution.x;
    else
        x0=cBest.x;
    end
    if iterm<0.50*StopL
        Candidate=selectCandidate(x0,Candidate_Num,varnum,iterm,StopL,LB,UB,pm,1);
    else
        Candidate=selectCandidate(x0,Candidate_Num,varnum,iterm,StopL,LB,UB,pm,2);
    end
    for i=1:Candidate_Num
         Candidate(i).fitness=fun(Candidate(i).x);
         f(i)=Candidate(i).fitness;
         if iterm==1
              Candidate(i).bTabu=0;
         end 
    end
    [a1,b1]=sort(f);
    for i=1:Candidate_Num
        for j=1:Tabu_Length
            a=norm(Candidate(i).x-TabuList(j,1:end-1));
            b=abs(Candidate(i).fitness-TabuList(j,end));
            if a<=0.1&&b<=0.005
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
    bestCandidate=Candidate(b1(1));   %最优候选解
    if bestCandidate.fitness<cBest.fitness
        m_solution=bestCandidate;
        for i=1:Tabu_Length-1
            TabuList(i,:)=TabuList(i+1,:);
        end
        TabuList(Tabu_Length,1:end-1)=bestCandidate.x;
        TabuList(Tabu_Length,end)=bestCandidate.fitness;
    else
        if balltabu==0
            for k=1:Candidate_Num
                if Candidate(k).bTabu==0
                    m_solution=Candidate(k);
                    for j=1:Tabu_Length-1
                       TabuList(j,:)=TabuList(j+1,:);
                    end
                    TabuList(Tabu_Length,1:end-1)=bestCandidate.x;
                    TabuList(Tabu_Length,end)=bestCandidate.fitness;
                    break
               end
            end
        else
            m_solution=bestCandidate;
            for j=1:Tabu_Length-1
                 TabuList(j,:)=TabuList(j+1,:);
            end
            TabuList(Tabu_Length,1:end-1)=bestCandidate.x;
            TabuList(Tabu_Length,end)=bestCandidate.fitness;
        end
    end
    if m_solution.fitness<cBest.fitness
        cBest=m_solution;
        cBest.index=iterm;
    end
end
