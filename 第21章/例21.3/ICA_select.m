function m_antibody=ICA_select(m_antibody)
antibodynum=size(m_antibody,2);
a=0.7;
P=zeros(1,antibodynum);
for i=1:antibodynum
   P(i)=a*m_antibody(i).pf+(1-a)*m_antibody(i).pd;
end
c=zeros(1,antibodynum);
for i=1:antibodynum
    if i==1
        c(i)=P(i);
    else
        c(i)=c(i-1)+P(i);
    end
end
c=c/c(antibodynum);
for i=1:antibodynum
    p=rand;
    index=1;
    while c(index)<p
        index=index+1;
    end
    newantibody(i)=m_antibody(index);
end
m_antibody=newantibody;