function y=lack_iter(A,x)   %��A��ȣ�xȱ�ٵ���Ŀ
x=unique(x);
m=length(x);
y1=[];
for i=1:m
    a=find(A==x(i));
    if ~isempty(a)
       y1=[y1 a];
    end
end
y=redu(A,y1,'c');