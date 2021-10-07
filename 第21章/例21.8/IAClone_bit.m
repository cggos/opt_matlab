function [best_x,fval]=IAClone_bit(fun,antibodynum,iterm_max,pm,nvar)        %�����ƿ�¡ѡ���㷨���󼫴�
for i=1:antibodynum
    m_antibody(i).x=rand(1,nvar)<0.5;
    m_antibody(i).density=0;
    m_antibody(i).pf=0;
    m_antibody(i).pd=0;
    m_antibody(i).fitness=fun(m_antibody(i).x);
    y(i)=m_antibody(i).fitness;
end
[a,b]=sort(y,'descend');      %�󼫴�
m_antibody=m_antibody(b);
antibody_bacterin=m_antibody(1:2);%�����
cbest=m_antibody(1);
for iterm=1:iterm_max
    cloneantibody=clone(m_antibody);%��¡
    cloneantibody=clonemutation_bit(cloneantibody,pm);%����
    for i=1:size(cloneantibody,2)
        cloneantibody(i).fitness=fun(cloneantibody(i).x);
        y1(i)=cloneantibody(i).fitness;
    end
    [a1,b1]=sort(y1,'descend');
    cloneantibody=cloneantibody(b1);
    m_antibody=cloneselection(cloneantibody,m_antibody);%ѡ����Դ��濹��Ŀ�¡��   
    antibody_bacterin=ICA_bacterin(m_antibody,antibody_bacterin);
    m_antibody=ICA_operator_bit(fun,m_antibody,antibody_bacterin);
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
