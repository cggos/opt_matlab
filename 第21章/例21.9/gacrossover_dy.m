function m_antibody=gacrossover_dy(m_antibody,px,xtype)%原来的基因是popsize×ndim矩阵
%ndim为基因的维数，popsize为种群规模
%交换算子
antibodynum=size(m_antibody,2);
NC=length(m_antibody(1).pattern);
halfpop=antibodynum/2;
if xtype==1%单点交换
    randlist=rand(halfpop,1);
    for i=1:halfpop
        a=(i*2)-1;
        xpo=a+1;
        if (randlist(i)<px)
            xpoint=ceil(rand*NC);
            temp=m_antibody(a).pattern(1:xpoint);
            m_antibody(a).pattern(1:xpoint)=m_antibody(xpo).pattern(1:xpoint);
            m_antibody(xpo).pattern(1:xpoint)=temp;
        end
    end
end
if xtype==2%两点交换
    randlist=rand(halfpop,1);
    for i=1:halfpop
        a=(i*2)-1;
        xpo=a+1;
        if (randlist(i)<px)
            xpoint=sort(ceil(rand(1,2).*NC));
            while xpoint(1)==xpoint(2)
                 xpoint=sort(ceil(rand(1,2).*NC));
            end
            temp=m_antibody(a).pattern(xpoint(1):xpoint(2));
            m_antibody(a).pattern(xpoint(1):xpoint(2))=m_antibody(xpo).pattern(xpoint(1):xpoint(2));
            m_antibody(xpo).pattern(xpoint(1):xpoint(2))=temp;
        end
    end
end
if xtype==3%多点交换
      for i=1:halfpop
        a=(i*2)-1;
        xpo=a+1;
         for j=1:NC
            test=rand;
            if test<px
               temp=m_antibody(a).pattern(j);
               m_antibody(a).pattern(j)=m_antibody(xpo).pattern(j);
               m_antibody(xpo).pattern(j)=temp;
            end
        end
     end
end