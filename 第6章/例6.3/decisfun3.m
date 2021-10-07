function u=decisfun3(k,x)
ww=[3,4,2,5];
if k==1 %物品1最多装载4件
   u=0:4;
   while x-ww(k)*u<0
   end
elseif k==2 %物品2最多装载3件
   u = 0:3;
   while x-ww(k)*u<0
   end
elseif k==3 %物品3最多装载6件
    u=0:6;
    while x-ww(k)*u<0
    end
else
    u=floor(x/ww(k));
end
u=u(:);