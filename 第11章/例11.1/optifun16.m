function y=optifun16(x)
g1=gcon1(x);
g2=gcon2(x);
punishment1=100;
punishment2=200;
punishment3=1000;
punishment4=3000;
q1=max(0,g1);
%q1=min(0,g1);
gamma1=1;
if q1<=1e-9
     theta1=0;
elseif q1>1e-9&&q1<0.001
     theta1=punishment1;
elseif q1>=0.001&&q1<0.1
     theta1=punishment2;
elseif q1>=0.1&&q1<1
     theta1=punishment3;
else 
     theta1=punishment4;
     gamma1=2;
end
lamda1=theta1*q1^gamma1;
q2=max(0,g2);
%q2=min(0,g2);
gamma2=1;
if q2<=1e-9
     theta2=0;
elseif q2>1e-9&&q2<0.001
     theta2=punishment1;
elseif q2>=0.001&&q2<0.1
     theta2=punishment2;
elseif q2>=0.1&&q2<1
     theta2=punishment3;
else 
     theta2=punishment4;
     gamma2=2;
end
lamda2=theta2*q2^gamma2;
y=100*(x(2)-x(1)^2)^2+(1-x(1))^2+(lamda1+lamda2);
   
function result=gcon1(x)
result=-x(1)-x(2)^2;
    
function result=gcon2(x)
result=-x(1)^2-x(2);  