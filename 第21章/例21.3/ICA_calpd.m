function m_antibody=ICA_calpd(m_antibody) %免疫中浓度概率计算
antibodynum=size(m_antibody,2);
m_antibody=ICA_densitysort(m_antibody);
T=(m_antibody(1).density+m_antibody(antibodynum).density)/2;
n=0;
pd1=0;
pd2=0;
 for i=1:antibodynum
    if m_antibody(i).density<T
       pd1=(1-(i-1)/antibodynum)/antibodynum;
       pd2=(1+((i-1)/antibodynum)*((i-1)/(antibodynum-(i-1))))/antibodynum;
       n=i;     %低浓度与高浓度的分界点
       break
    else
        continue
    end
 end
 for i=1:antibodynum
     if i<n
         m_antibody(i).pd=pd1;
     else
         m_antibody(i).pd=pd2;
     end
 end
         