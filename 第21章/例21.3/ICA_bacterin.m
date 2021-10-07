function antibody_bacterin=ICA_bacterin(m_antibody,antibody_bacterin)    %免疫算法，构建疫苗
antibodynum=size(m_antibody,2);
n=size(antibody_bacterin,2);
if n~=ceil(0.2*antibodynum)%再加入两个疫苗
     for i=1:n+2
         if i<=n
             temp1(i)=antibody_bacterin(i);
         else
             temp1(i)=m_antibody(i-n);
         end
     end
     antibody_bacterin=temp1;
     for i=1:n+1   %使疫苗按从大到小排列
        for j=i:n+2
           if antibody_bacterin(i).fitness<antibody_bacterin(j).fitness
                 temp=antibody_bacterin(i);
                 antibody_bacterin(i)=antibody_bacterin(j);
                 antibody_bacterin(j)=temp;
           end
        end
     end       
else   %替换
   temp(1)=m_antibody(1);temp(2)=m_antibody(2);             
   for i=1:2
       for j=1:n
           if temp(i).fitness>antibody_bacterin(j).fitness
              antibody_bacterin(j)=temp(i);
              break
           end
       end
   end
end   