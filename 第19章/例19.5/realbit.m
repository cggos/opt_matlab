function b=realbit(x)   %实数转化成二进制
n=length(x);
b=zeros(1,n);
for i=1:n
    if 1/(1+exp(-x(i)))<=0.5
        b(i)=0;
    else
        b(i)=1;
    end
end