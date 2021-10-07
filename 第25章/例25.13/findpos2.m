function y=findpos2(x,x1)   %x1各元素在x中的位置
num=length(x1);
for i=1:num
    a=find(x==x1(i));
    if isempty(a)
        y(i)=[];
    else
        y(i)=a;
    end
end