function [best_x,fval]=IAGA(fun,antibodynum,iterm_max,pc,pm,LB,UB)       %√‚“ﬂ“≈¥´À„∑®
numvar=size(LB,1);
for i=1:antibodynum
    m_antibody(i).x=LB'+rand(1,numvar).*(UB-LB)';    
    m_antibody(i).density=0;
    m_antibody(i).pf=0;
    m_antibody(i).pd=0;
    m_antibody(i).fitness=fun(m_antibody(i).x);
    y(i)=m_antibody(i).fitness;
end
[a,b]=sort(y,'descend');    %«Ûº´¥Û
m_antibody=m_antibody(b);  %≈≈–Ú
antibody_bacterin=m_antibody(1:2);%“ﬂ√Áø‚
cbest=m_antibody(1);
for iterm=1:iterm_max
    m_antibody=gacrossover(m_antibody,pc,[],'r',LB,UB);
    m_antibody=gamutation(m_antibody,pm,LB,UB,iterm,iterm_max);   
    for i=1:antibodynum
         m_antibody(i).fitness=fun(m_antibody(i).x);
         y(i)=m_antibody(i).fitness;
    end
    [a,b]=sort(y,'descend');
    m_antibody=m_antibody(b);  %≈≈–Ú
    antibody_bacterin=ICA_bacterin(m_antibody,antibody_bacterin);
    m_antibody=ICA_operator(fun,m_antibody,antibody_bacterin,LB,UB);
    m_antibody=ICA_caldensity(m_antibody);
    m_antibody=ICA_calpf(m_antibody);
    m_antibody=ICA_calpd(m_antibody);
    m_antibody=ICA_select(m_antibody);
    for i=1:antibodynum
         m_antibody(i).fitness=fun(m_antibody(i).x);
         y(i)=m_antibody(i).fitness;
    end
    [a,b]=sort(y,'descend');
    if cbest.fitness<a(1)
        cbest.x=m_antibody(b(1)).x;
        cbest.fitness=a(1);
    end
end
best_x=cbest.x;
fval=cbest.fitness;





    
