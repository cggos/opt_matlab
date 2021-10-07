function cBest=TSOATSP(city,StopL,Candidate_Num,Tabu_Length)
%cityΪ�������꣬itermax����������,L���ɱ���
[m,Distance_List]=city2d(city);
%m=size(city,1);
PlaceNum=m;                   %���еص�����
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
TabuList=zeros(Tabu_Length,PlaceNum);   	%���ɱ�(tabu list)
m_solution.route=randperm(PlaceNum);      	%������ɳ�ʼ������;
m_solution.fitness=value(m_solution.route,Distance_List);
m_solution.idex=1;
cBest=m_solution;               %���Ž�
for iterm=1:StopL
    A=zeros(Candidate_Num,2);   
    A=ceil(PlaceNum*rand(Candidate_Num,2));  	    %�������Ԫ��ֵȡֵ��ΧΪ[0,PlaceNum]��2ά�����������ɺ�ѡ��
    for i=1:Candidate_Num             %��ѡ��
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
        Si(i,:)=S;                                  %���ú�ѡ���ÿһ��SΪ��ʼ�������
        Si(i,[A(i,1),A(i,2)])=S([A(i,2),A(i,1)]);   %������i�е�����λ�õ�����λ�õı����A����
        %�����i����ѡ����ܾ���ֵ
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
                Candidate(i).bTabu=1;     %�ڽ��ɱ���
                break
            else
                Candidate(i).bTabu=0;
            end
        end
    end
    balltabu=1;    %�ж��Ƿ��ڽ��ɱ���
    for i=1:Candidate_Num
       if Candidate(i).bTabu==0     %���ڽ��ɱ���
           balltabu=0;
           break
       end
    end
    bestCandidate=Candidate(1);   %���ź�ѡ��
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
Route=cBest.route;        %�����������Ž�Ļ���·��
TSPplot(city,Route);      %��ͼ
%gname;



