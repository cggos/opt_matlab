function m_antibody=ICA_operator_bit(fun,m_antibody,antibody_bacterin)%×¢ÈëÒßÃç
antibodynum=size(m_antibody,2);
n=size(antibody_bacterin,2);
cfitness_b=zeros(1,n);
cfitness_a=zeros(1,antibodynum);
px=0.6;
NC=length(m_antibody(1).x);
for i=1:n
    if i==1
        cfitness_b(i)=antibody_bacterin(i).fitness;
    else
        cfitness_b(i)=cfitness_b(i-1)+antibody_bacterin(i).fitness;
    end
end
cfitness_b=cfitness_b/cfitness_b(n);
for i=1:antibodynum
    if i==1
        cfitness_a(i)=m_antibody(i).fitness;
    else
        cfitness_a(i)=cfitness_a(i-1)+m_antibody(i).fitness;
    end
end
cfitness_a=cfitness_a/cfitness_a(antibodynum);
for i=1:n
    p=rand;
    if p<cfitness_b(i)
       index=i;    %ÕÒµ½ÒßÃç
       break;
    end
end
if index==n|isempty(index)
    index=1;
end
index1=[];
for i=1:antibodynum
    pick=rand;
    for j=1:antibodynum
        if pick<cfitness_a(j)
            index1=[index1 i];    %ÕÒµ½¿¹Ìå
            break;
        end
    end
end
m=length(index1);
for i=1:m
    temp=m_antibody(index1(i));    
    if rand<px
       xpoint=ceil(rand*NC);
       if xpoint>=NC/2
           temp.x(1:xpoint)=antibody_bacterin(index).x(1:xpoint);
       else
           temp.x(xpoint:end)=antibody_bacterin(index).x(xpoint:end);
       end
    end
    y=fun(temp.x);
    if y>m_antibody(index1(i)).fitness
        m_antibody(index1(i)).x=temp.x;
        m_antibody(index1(i)).fitness=y;
    end
end
m_antibody(end)=antibody_bacterin(index);
m_antibody(end-1)=antibody_bacterin(2);
for i=1:antibodynum-1
    for j=i:antibodynum
        if m_antibody(i).fitness<m_antibody(j).fitness
            temp=m_antibody(i);
            m_antibody(i)=m_antibody(j);
            m_antibody(j)=temp;
        end
    end
end