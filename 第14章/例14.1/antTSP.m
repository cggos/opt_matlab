function [Shortest_Route,Shortest_Length]=antTSP(city,iter_max,m,Alpha,Beta,Rho,Q)
%city：n个城市的坐标，n×2的矩阵
% iter_max ：最大迭代次数
% m： 蚂蚁个数
%Alpha ：表征信息素重要程度的参数
%Beta ：表征启发式因子重要程度的参数
% Rho ：信息素蒸发系数
% Q ：信息素增加强度系数
% R_best 各代最佳路线
[n,D]=city2d(city);
Eta=1./D;          %Eta为能见度因数，这里设为距离的倒数
Tau=ones(n,n);     %Tau为信息素矩阵
Tabu=zeros(m,n);   %存储并记录路径的生成
nC=1;             %计数器
R_best=zeros(iter_max,n);
L_best=inf.*ones(iter_max,1);%%各代最佳路线及长度
while nC<=iter_max    
  Tabu(:,1)=ceil(rand(m,1).*n);
   %m只蚂蚁按概率函数选择下一座城市，完成各自的周游
  for j=2:n
    for i=1:m
        visited=Tabu(i,1:(j-1));  %已访问的城市
        J=zeros(1,(n-j+1)); %待访问的城市
        P=J; Jc=1; %待访问城市的选择概率分布
        for k=1:n
           if length(find(visited==k))==0
               J(Jc)=k;Jc=Jc+1;
           end
       end
       for k=1:length(J)
           P(k)=(Tau(visited(end),J(k))^Alpha)*(Eta(visited(end),J(k))^Beta);%待选城市概率分布
       end
       P=P/(sum(P)); 
       Pcum=cumsum(P);
       Select=find(Pcum>=rand);%按概率选取下一个城市
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
  %记录本次迭代
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
%更新信息素，本例采用的蚁周系统　
  Delta_Tau=zeros(n,n);
  for i=1:m
      for j=1:(n-1)
          Delta_Tau(Tabu(i,j),Tabu(i,j+1))=Delta_Tau(Tabu(i,j),Tabu(i,j+1))+Q/L(i);
      end
      Delta_Tau(Tabu(i,n),Tabu(i,1))=Delta_Tau(Tabu(i,n),Tabu(i,1))+Q/L(i);
  end
  Tau=(1-Rho).*Tau+Delta_Tau;
  %禁忌表清零
  Tabu=zeros(m,n);
end
%输出结果
Pos=find(L_best==min(L_best));
Shortest_Route=R_best(Pos(1),:);
Shortest_Length=L_best(Pos(1));
TSPplot(city,Shortest_Route);