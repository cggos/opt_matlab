function f=subfun5(k,x1,x2,u1,u2)
A=[0 1 3 6;4 5 6 7;5 6 7 8;6 7 8 9];
B=[0 2 4 6;1 4 6 7;4 6 8 9;6 8 10 11];
C=[0 3 5 8;2 5 7 9;4 7 9 11;6 9 11 13];
if k==1
    f=-A(u1+1,u2+1);%求最大值转化为求
elseif k==2
    f=-B(u1+1,u2+1);
else
    f=-C(u1+1,u2+1);
end


    
