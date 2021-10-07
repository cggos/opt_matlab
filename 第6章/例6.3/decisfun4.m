function u=decisfun4(k,x)
ww=[3,4,2,5];
if k==1 %物
    u=1:4;
    while x-ww(k)*u<0
    end
elseif k==2 %物品2至少装载1件，最多装载3件
    u=1:3;
    while x-ww(k)*u<0
    end
elseif k==3 %物品3最多装载6件
    u=0:6;
    while x-ww(k)*u < 0
    end
else
    u=floor(x/ww(k)); %所装载物品的数量为该阶段的最大承载量除以物品单位重量
end
u=u(:);
