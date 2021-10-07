function m_antibody=ICA_caldensity(m_antibody)  %免疫算法计算抗体浓度
antibodynum=size(m_antibody,2);
m_antibody=ICA_fitsort(m_antibody);
b1=m_antibody(1);
n1=1;
for i=1:antibodynum
    %if (m_antibody(i).fitness/b1.fitness)<1.02
    if abs(m_antibody(i).fitness-b1.fitness)<0.5
        continue
    else
        b1=m_antibody(i);
        n2=i;
        num=n2-n1;
        for j=1:num
            m_antibody(n1+j-1).density=num/antibodynum;
        end
        n1=i;
    end
    if i==antibodynum
        n2=i+1;
        num=n2-n1;
        for j=1:num
            m_antibody(n1+j-1).density=num/antibodynum;
        end
    end
end

      
        
    