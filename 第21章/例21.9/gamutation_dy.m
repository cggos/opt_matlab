function m_antibody=gamutation_dy(m_antibody,pm)
%±äÒìËã×Ó
antibodynum=size(m_antibody,2);
NC=length(m_antibody(1).pattern);
classnum=length(unique(m_antibody(1).pattern));
for i=1:antibodynum
    for j=1:NC
        test=rand;
        if test<pm
            b=m_antibody(i).pattern(j);
            a=ceil(classnum*rand);
            while a==b
               a=ceil(classnum*rand);
            end
             m_antibody(i).pattern(j)=a;
        end
    end
end
            