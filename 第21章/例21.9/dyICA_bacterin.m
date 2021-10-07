function antibody_bacterin=dyICA_bacterin(m_antibody,antibody_bacterin,bacterinnum)  %∂ØÃ¨∏¸ªª“ﬂ√Á 
for i=1:bacterinnum
   for j=1:bacterinnum
     if m_antibody(i).fitness>antibody_bacterin(j).fitness
         antibody_bacterin(j).fitness=m_antibody(i).fitness;
     end
   end
end   
for i=1:bacterinnum-1    %≈≈–Ú
    for j=i:bacterinnum
       if antibody_bacterin(i).fitness<antibody_bacterin(j).fitness
            temp=antibody_bacterin(i).fitness;
            antibody_bacterin(i).fitness=antibody_bacterin(j).fitness;
            antibody_bacterin(j).fitness=temp;
        end
    end
end
