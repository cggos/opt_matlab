function cloneantibody=clone(m_antibody)  %ПЫТЁЫузг
antibodynum=size(m_antibody,2);
cloneantibodynum=0;
for i=1:antibodynum
    n=floor(sqrt(antibodynum/i));
    cloneantibodynum=cloneantibodynum+n;
    for j=1:n
        cloneantibody(cloneantibodynum-n+j)=m_antibody(i);
    end
end
