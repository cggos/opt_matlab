function y=center_dis(m_center)%«Û¿‡º‰æ‡¿Î
[centern,c]=size(m_center);
y=0;   
for i=1:centern
    for j=i+1:centern
        for k=1:c
          y=y+(m_center(i,k)-m_center(j,k))^2;
        end
    end
end

