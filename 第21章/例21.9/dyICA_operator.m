function m_antibody=dyICA_operator(fun,m_antibody,antibody_bacterin,bacterinnum)
antibodynum=size(m_antibody,2);
NC=length(m_antibody(1).pattern);
cfitness_a=zeros(1,antibodynum);
for i=1:antibodynum
    if i==1
        cfitness_a(i)=m_antibody(i).fitness;
    else
        cfitness_a(i)=cfitness_a(i-1)+m_antibody(i).fitness;
    end
end
cfitness_a=cfitness_a/cfitness_a(antibodynum);
for i=1:antibodynum
    pick=rand;
    index=1;
    while cfitness_a(index)<pick
        index=index+1;
    end
    newantibody=m_antibody(index);  %选择抗体
    xpoint=zeros(1,2);
    while xpoint(1)==xpoint(2)
        xpoint=sort(ceil(NC*rand(1,2)));
    end
    n=ceil(bacterinnum*rand);%随机选择疫苗
    newantibody.pattern(xpoint(1):xpoint(2))=antibody_bacterin(n).pattern(xpoint(1):xpoint(2));
    newantibody.fitness=fun(newantibody.pattern);
    if newantibody.fitness>m_antibody(i).fitness
        m_antibody(i)=newantibody;
    end
end
for i=1:antibodynum-1
    for j=i:antibodynum
        if m_antibody(i).fitness<m_antibody(j).fitness
            temp=m_antibody(i);
            m_antibody(i)=m_antibody(j);
            m_antibody(j)=temp;
        end
    end
end
