function y=lack_iter(A,x)   %与A相比，x缺少的项目
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