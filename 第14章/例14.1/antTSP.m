function [Shortest_Route,Shortest_Length]=antTSP(city,iter_max,m,Alpha,Beta,Rho,Q)
%city��n�����е����꣬n��2�ľ���
% iter_max ������������
% m�� ���ϸ���
%Alpha ��������Ϣ����Ҫ�̶ȵĲ���
%Beta ����������ʽ������Ҫ�̶ȵĲ���
% Rho ����Ϣ������ϵ��
% Q ����Ϣ������ǿ��ϵ��
% R_best �������·��
[n,D]=city2d(city);
Eta=1./D;          %EtaΪ�ܼ���������������Ϊ����ĵ���
Tau=ones(n,n);     %TauΪ��Ϣ�ؾ���
Tabu=zeros(m,n);   %�洢����¼·��������
nC=1;             %������
R_best=zeros(iter_max,n);
L_best=inf.*ones(iter_max,1);%%�������·�߼�����
while nC<=iter_max    
  Tabu(:,1)=ceil(rand(m,1).*n);
   %mֻ���ϰ����ʺ���ѡ����һ�����У���ɸ��Ե�����
  for j=2:n
    for i=1:m
        visited=Tabu(i,1:(j-1));  %�ѷ��ʵĳ���
        J=zeros(1,(n-j+1)); %�����ʵĳ���
        P=J; Jc=1; %�����ʳ��е�ѡ����ʷֲ�
        for k=1:n
           if length(find(visited==k))==0
               J(Jc)=k;Jc=Jc+1;
           end
       end
       for k=1:length(J)
           P(k)=(Tau(visited(end),J(k))^Alpha)*(Eta(visited(end),J(k))^Beta);%��ѡ���и��ʷֲ�
       end
       P=P/(sum(P)); 
       Pcum=cumsum(P);
       Select=find(Pcum>=rand);%������ѡȡ��һ������
       if  isempty(Select)
           Tabu(i,j)=round(1+(n-1)*rand);
       else
           next_visit=J(Select(1));Tabu(i,j)=next_visit;
       end
    end
  end
  if nC>=2
      Tabu(1,:)=R_best(nC-1,:);
  end
  %��¼���ε���
  L=zeros(m,1);
  for i=1:m
      R=Tabu(i,:);
      for j=1:(n-1);
         L(i)=L(i)+D(R(j),R(j+1));
      end 
      L(i)=L(i)+D(R(1),R(n));
  end 
  L_best(nC)=min(L);
  pos=find(L==L_best(nC));
  R_best(nC,:)=Tabu(pos(1),:);
  nC=nC+1;
%������Ϣ�أ��������õ�����ϵͳ��
  Delta_Tau=zeros(n,n);
  for i=1:m
      for j=1:(n-1)
          Delta_Tau(Tabu(i,j),Tabu(i,j+1))=Delta_Tau(Tabu(i,j),Tabu(i,j+1))+Q/L(i);
      end
      Delta_Tau(Tabu(i,n),Tabu(i,1))=Delta_Tau(Tabu(i,n),Tabu(i,1))+Q/L(i);
  end
  Tau=(1-Rho).*Tau+Delta_Tau;
  %���ɱ�����
  Tabu=zeros(m,n);
end
%������
Pos=find(L_best==min(L_best));
Shortest_Route=R_best(Pos(1),:);
Shortest_Length=L_best(Pos(1));
TSPplot(city,Shortest_Route);