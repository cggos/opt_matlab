function [best_route,fval]=ICA_TSP(city,antibodynum,iterm_max,pc,pm)%√‚“ﬂÀ„∑®
[NC,d]=city2d(city);
for i=1:antibodynum
    m_antibody(i).route=randperm(NC);
    m_antibody(i).density=0;
    m_antibody(i).pf=0;
    m_antibody(i).pd=0;
    m_antibody(i).fitness=1/value(m_antibody(i).route,d);
    y(i)=m_antibody(i).fitness;
end
[a,b]=sort(y,'descend');
m_antibody=m_antibody(b);%“—æ≠≈≈–Ú
antibody_bacterin=m_antibody(1:2);%“ﬂ√Áø‚
cbest=m_antibody(1);
for iterm=1:iterm_max
    m_antibody=gacrossover_TSP(m_antibody,pc);%¡Ωµ„Ωª≤Ê
    m_antibody=gamutation_TSP(m_antibody,pm);
    for i=1:antibodynum
        m_antibody(i).fitness=1/value(m_antibody(i).route,d);
        y(i)=m_antibody(i).fitness;
    end
    [a,b]=sort(y,'descend');
    m_antibody=m_antibody(b);  %≈≈–Ú
    antibody_bacterin=ICA_bacterin(m_antibody,antibody_bacterin);
    m_antibody=ICA_operator_TSP(m_antibody,antibody_bacterin,d);
    m_antibody=ICA_caldensity(m_antibody);
    m_antibody=ICA_calpf(m_antibody);
    m_antibody=ICA_calpd(m_antibody);
    m_antibody=ICA_select(m_antibody);
    for i=1:antibodynum
        m_antibody(i).fitness=1/value(m_antibody(i).route,d);
        y(i)=m_antibody(i).fitness;
    end
    [a,b]=sort(y,'descend');
    if cbest.fitness<a(1)
        cbest.route=m_antibody(b(1)).route;
        cbest.fitness=a(1);
    end
end
best_route=cbest.route;
fval=1/cbest.fitness;
TSPplot(city,best_route);



    
