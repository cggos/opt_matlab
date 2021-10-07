function m_antibody=gacrossover_TSP(m_antibody,px)%原来的基因是popsize×ndim矩阵
%ndim为基因的维数，popsize为种群规模
%交换算子
antibodynum=size(m_antibody,2);
NC=length(m_antibody(1).route);
halfpop=antibodynum/2;
randlist=rand(halfpop,1);
for k=1:halfpop  
   a=(k*2)-1;
   xpo=a+1;
   if (randlist(k)<px)
     xpoint=sort(ceil(rand(1,2).*NC));
     while xpoint(2)==xpoint(1)
         xpoint=sort(ceil(rand(1,2).*NC));
     end
     if xpoint(2)==NC
          temp1=my_delrow(m_antibody(a).route,m_antibody(xpo).route(xpoint(1):xpoint(2)));
          temp2=my_delrow(m_antibody(xpo).route,m_antibody(a).route(xpoint(1):xpoint(2)));
          m_antibody(a).route(1:xpoint(1)-1)=temp2;
          m_antibody(xpo).route(1:xpoint(1)-1)=temp1;
     elseif xpoint(1)==1
          temp1=my_delrow(m_antibody(a).route,m_antibody(xpo).route(xpoint(1):xpoint(2)));
          temp2=my_delrow(m_antibody(xpo).route,m_antibody(a).route(xpoint(1):xpoint(2)));
          m_antibody(a).route(xpoint(2)+1:NC)=temp2;
          m_antibody(xpo).route(xpoint(2)+1:NC)=temp1;
     elseif xpoint(1)==1&&xpoint(2)==NC
         continue
     else         
          m1=xpoint(1)-1;
          temp1=my_delrow(m_antibody(a).route,m_antibody(xpo).route(xpoint(1):xpoint(2)));
          temp2=my_delrow(m_antibody(xpo).route,m_antibody(a).route(xpoint(1):xpoint(2)));
          m_antibody(a).route(1:xpoint(1)-1)=temp2(1:m1);
          m_antibody(a).route(xpoint(2)+1:NC)=temp2(m1+1:end);
          m_antibody(xpo).route(1:xpoint(1)-1)=temp1(1:m1);
          m_antibody(xpo).route(xpoint(2)+1:NC)=temp1(m1+1:end);
     end
   end
end



