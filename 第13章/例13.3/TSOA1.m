function cBest=TSOA1(fun)   %������ѡ��
prompt={'��������';'��ѡ�����';'���ɱ���';'����ά��'};
name='�����㷨������';
defaultanswer={'1000','40','10','[]'};
answer=inputdlg(prompt,name,1,defaultanswer);
StopL=str2num(answer{1});
Candidate_Num=str2num(answer{2});
Tabu_Length=str2num(answer{3});
m=str2num(answer{4});
TabuList=zeros(Tabu_Length,m);   	%���ɱ�(tabu list)
m_solution.pattern=rand(1,m)<0.5;      	%������ɳ�ʼ������;
%m_solution.fitness=fun(m_solution.pattern,data,target);
m_solution.fitness=fun(m_solution.pattern);
m_solution.idex=1;
cBest=m_solution;               %���Ž�
for iterm=1:StopL
    A=zeros(Candidate_Num,2);   
    A=ceil(m*rand(Candidate_Num,2));  	    %�������Ԫ��ֵȡֵ��ΧΪ[0,m]��2ά�����������ɺ�ѡ��
    for i=1:Candidate_Num    %��ѡ��
       while isequal(A(i,1),A(i,2))     
           A(i,:)=ceil(m*rand(1,2));  
       end       
       for j=1:Candidate_Num
          if j~=i
             while isequal(A(i,:),A(j,:))     
                A(j,:)=ceil(m*rand(1,2));  
             end
          end
       end
     end
    S=m_solution.pattern;
    Si=zeros(Candidate_Num,m);
    for i=1:Candidate_Num
        Si(i,:)=S;                                  %���ú�ѡ���ÿһ��SΪ��ʼ�������
        Si(i,[A(i,1),A(i,2)])=S([A(i,2),A(i,1)]);   %������i�е�����λ�õ�����λ�õı����A����
        %�����i����ѡ����ܾ���ֵ
        Total_D(i)=fun(Si(i,:));
    end
    [a,b]=sort(Total_D);
    for i=1:Candidate_Num
        Candidate(i).fitness=Total_D(b(i));
        Candidate(i).pattern=Si(b(i),:);
        if iterm==1
            Candidate(i).bTabu=0;
        end 
    end
    for i=1:Candidate_Num
        for j=1:Tabu_Length
            if isequal(Candidate(i).pattern,TabuList(j,:))
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
        TabuList(Tabu_Length,:)=bestCandidate.pattern;
    else
        if balltabu==0
            for k=1:Candidate_Num
                if Candidate(k).bTabu==0
                    m_solution=Candidate(k);
                    for j=1:Tabu_Length-1
                       TabuList(j,:)=TabuList(j+1,:);
                    end
                    TabuList(Tabu_Length,:)=bestCandidate.pattern;
                    break
               end
            end
        else
            m_solution=bestCandidate;
            for j=1:Tabu_Length-1
                 TabuList(j,:)=TabuList(j+1,:);
            end
            TabuList(Tabu_Length,:)=bestCandidate.pattern;
        end
    end
    if m_solution.fitness<cBest.fitness
        cBest=m_solution;
        cBest.index=iterm;
    end
end




