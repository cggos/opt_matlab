function u=decisfun3(k,x)
ww=[3,4,2,5];
if k==1 %��Ʒ1���װ��4��
   u=0:4;
   while x-ww(k)*u<0
   end
elseif k==2 %��Ʒ2���װ��3��
   u = 0:3;
   while x-ww(k)*u<0
   end
elseif k==3 %��Ʒ3���װ��6��
    u=0:6;
    while x-ww(k)*u<0
    end
else
    u=floor(x/ww(k));
end
u=u(:);