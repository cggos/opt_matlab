function u=decisfun4(k,x)
ww=[3,4,2,5];
if k==1 %��
    u=1:4;
    while x-ww(k)*u<0
    end
elseif k==2 %��Ʒ2����װ��1�������װ��3��
    u=1:3;
    while x-ww(k)*u<0
    end
elseif k==3 %��Ʒ3���װ��6��
    u=0:6;
    while x-ww(k)*u < 0
    end
else
    u=floor(x/ww(k)); %��װ����Ʒ������Ϊ�ý׶ε���������������Ʒ��λ����
end
u=u(:);
