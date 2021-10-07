function [best_pattern,fval]=IA_cluster(fun,antibodynum,iterm_max,pc,pm,classnum,NC)     %动态疫苗免疫算法,聚类
for i=1:antibodynum
    m_antibody(i).pattern=unidrnd(classnum,[1 NC]);    
    m_antibody(i).density=0;
    m_antibody(i).pf=0;
    m_antibody(i).pd=0;
    m_antibody(i).fitness=fun(m_antibody(i).pattern);
    y(i)=m_antibody(i).fitness;
end
[a,b]=sort(y,'descend');
m_antibody=m_antibody(b);%已经排序
bacterinnum=ceil(antibodynum*0.05);
antibody_bacterin=m_antibody(1:bacterinnum);%疫苗库，数量固定为3个
cbest=m_antibody(1);
for iterm=1:iterm_max
    m_antibody=gacrossover_dy(m_antibody,pc,1);%两点交叉
    m_antibody=gamutation_dy(m_antibody,pm);
    for i=1:antibodynum
        m_antibody(i).fitness=fun(m_antibody(i).pattern);
        y(i)=m_antibody(i).fitness;
    end
    [a,b]=sort(y,'descend');
    m_antibody=m_antibody(b);%已经排序
    antibody_bacterin=dyICA_bacterin(m_antibody,antibody_bacterin,bacterinnum);
    m_antibody=dyICA_operator(fun,m_antibody,antibody_bacterin,bacterinnum);
    m_antibody=ICA_caldensity(m_antibody);
    m_antibody=ICA_calpf(m_antibody);
    m_antibody=ICA_calpd(m_antibody);
    m_antibody=ICA_select(m_antibody);
    for i=1:antibodynum
         m_antibody(i).fitness=fun(m_antibody(i).pattern);
         y(i)=m_antibody(i).fitness;
    end
    [a,b]=sort(y,'descend');
    if cbest.fitness<a(1)
       cbest.x=m_antibody(b(1)).pattern;
       cbest.fitness=a(1);
    end
end
best_pattern=cbest.pattern;
fval=1/cbest.fitness;

