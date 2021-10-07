function m_antibody=ICA_calpf(m_antibody)  %«Û≈®∂»∏≈¬ 
antibodynum=size(m_antibody,2);
totalfitness=0;
for i=1:antibodynum
    totalfitness=totalfitness+m_antibody(i).fitness;
end
for i=1:antibodynum
    m_antibody(i).pf=m_antibody(i).fitness/totalfitness;
end