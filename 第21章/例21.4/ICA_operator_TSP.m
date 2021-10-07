function m_antibody=ICA_operator_TSP(m_antibody,antibody_bacterin,d)%×¢ÈëÒßÃç
antibodynum=size(m_antibody,2);
NC=length(m_antibody(1).route);
n=size(antibody_bacterin,2);
cfitness_b=zeros(1,n);
cfitness_a=zeros(1,antibodynum);
px=0.6;
for i=1:n
    if i==1
        cfitness_b(i)=antibody_bacterin(i).fitness;
    else
        cfitness_b(i)=cfitness_b(i-1)+antibody_bacterin(i).fitness;
    end
end
cfitness_b=cfitness_b/cfitness_b(n);
for i=1:antibodynum
    if i==1
        cfitness_a(i)=m_antibody(i).fitness;
    else
        cfitness_a(i)=cfitness_a(i-1)+m_antibody(i).fitness;
    end
end
cfitness_a=cfitness_a/cfitness_a(antibodynum);
for i=1:n
    p=rand;
    p=p-cfitness_b(i);
       if p<0
          index=i;    %ÕÒµ½ÒßÃç
          break;
       end
end
if index==n|isempty(index)
    index=1;
end
index1=[];
for i=1:antibodynum
    pick=rand;
    for j=1:antibodynum
        pick=pick-cfitness_a(j);
        if pick<0
            index1=[index1 i];    %ÕÒµ½¿¹Ìå
            break;
        end
    end
end
m=length(index1);
for i=1:m
   test=rand;
      if test<px
          xpoint=ceil(rand*NC);
          temp=m_antibody(index1(i));
          if xpoint>ceil(NC/2)
              temp.route(1:xpoint)=antibody_bacterin(index).route(1:xpoint);
              temp.route=isin_TSP(temp.route);
          else
              temp.route(xpoint:end)=antibody_bacterin(index).route(xpoint:end);
              temp.route=isin_TSP(temp.route);
          end
          y=1/value(temp.route,d);        
          if y>temp.fitness
             m_antibody(index1(i))=temp;
          end        
      end  
end
m_antibody(end)=antibody_bacterin(index);
m_antibody(end-1)=antibody_bacterin(2);
for i=1:antibodynum-1
   for j=i:antibodynum
      if m_antibody(i).fitness<m_antibody(j).fitness
         temp1=m_antibody(i);
         m_antibody(i)=m_antibody(j);
         m_antibody(j)=temp1;
      end
   end
end
  
  
  
   

