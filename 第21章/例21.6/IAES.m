function [best_x,fval]=IAES(fun,antibodynum,iterm_max,LB,UB)       %免疫规划算法,求极大
nvar=size(LB,1);   %变量的个数
sigma=zeros(antibodynum,nvar);
lamda=7*antibodynum;
for i=1:antibodynum
    m_antibody(i).x=LB'+rand(1,nvar).*(UB-LB)';   
    m_antibody(i).density=0;
    m_antibody(i).pf=0;
    m_antibody(i).pd=0;
    m_antibody(i).fitness=fun(m_antibody(i).x);
    y(i)=m_antibody(i).fitness;   
end
for j=1:nvar
   sigma(:,j)=3.0.*ones(antibodynum,1);
end
[a,b]=sort(y,'descend');    %求极大
m_antibody=m_antibody(b);  %排序
antibody_bacterin=m_antibody(1:2);%疫苗库
cbest=m_antibody(1);
for iterm=1:iterm_max
    [newchome,newsigma]=IAES_recombination1(m_antibody,sigma,lamda,3);
    [newchome,newsigma]=IAES_mutation(newchome,newsigma,LB,UB);
    for i=1:antibodynum
        newchome(i).fitness=fun(newchome(i).x);
    end
    [m_antibody,sigma]=IAES_select1(m_antibody,newchome,sigma,newsigma,2);
    for j=1:antibodynum
        m_antibody(j).fitness=fun(m_antibody(j).x);
        y(j)=m_antibody(j).fitness; 
    end   
    [a,b]=sort(y,'descend');
    m_antibody=m_antibody(b);  %排序
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





    
