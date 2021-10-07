function y=optifun15_1(x,LB,UB,CodeLen)
temp=Dec(LB,UB,x,CodeLen);
y=sin(temp(1))*sin(temp(2))/temp(1)/temp(2);