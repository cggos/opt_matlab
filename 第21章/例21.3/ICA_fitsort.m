function m_antibody=ICA_fitsort(m_antibody)  %�����㷨��Ũ�ȵ����򣨴Ӵ�С��
antibodynum=size(m_antibody,2);
for i=1:antibodynum-1
    for j=i:antibodynum
        if m_antibody(i).fitness<m_antibody(j).fitness
            temp=m_antibody(i);
            m_antibody(i)=m_antibody(j);
            m_antibody(j)=temp;
        end
    end
end