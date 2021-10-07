function [best_x,fval]=IAEP(fun,antibodynum,iterm_max,LB,UB,type)       %免疫规划算法,求极大
nvar=size(LB,1);   %变量的个数
for i=1:antibodynum
    m_antibody(i).x=LB'+rand(1,nvar).*(UB-LB)';   
    m_antibody(i).density=0;
    m_antibody(i).pf=0;
    m_antibody(i).pd=0;
    m_antibody(i).fitness=fun(m_antibody(i).x);
    m_antibody(i).sigma=sqrt(m_antibody(i).fitness);
    y(i)=m_antibody(i).fitness;    
end
for i=1:nvar
    x=[];
    for j=1:antibodynum
        x=[x;m_antibody(j).x(i)];
    end
    m_antibody_Var(i)=var(x);  
end
[a,b]=sort(y,'descend');    %求极大
m_antibody=m_antibody(b);  %排序
antibody_bacterin=m_antibody(1:2);%疫苗库
cbest=m_antibody(1);
for iterm=1:iterm_max
    m_antibody=gaussmution(m_antibody,m_antibody_Var,LB,UB,type);
    for i=1:antibodynum
         m_antibody(i).fitness=fun(m_antibody(i).x);
         y(i)=m_antibody(i).fitness;
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
         m_antibody(i).sigma=sqrt(m_antibody(i).fitness);
         y(i)=m_antibody(i).fitness;
    end
    [a,b]=sort(y,'descend');
    if cbest.fitness<a(1)
        cbest.x=m_antibody(b(1)).x;
        cbest.fitness=a(1);
    end
    for i=1:nvar
       x=[];
       for j=1:antibodynum
          x=[x;m_antibody(j).x(i)];
       end
    m_antibody_Var(i)=var(x);  
    end
end
best_x=cbest.x;
fval=cbest.fitness;





    
