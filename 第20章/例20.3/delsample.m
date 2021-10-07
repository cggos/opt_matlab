function y=delsample(x)   %去掉重复样本
r=size(x,1);
if r==1
    y=x;
else
    y=[];
    while 1
        temp=x(1,:);
        y=[y;temp]; 
        s=[];
        for i=1:r
           if (x(i,:)-temp)<1e-4
              s=[s i];
           end
        end
        x=redu(x,s,'r');
        if isempty(x)
           break
        end
        r=size(x,1);
    end
end