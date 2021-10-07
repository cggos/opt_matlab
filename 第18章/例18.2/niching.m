function [new,new1,k]=niching(pop,LB,UB,n)   %生成小生境 ,pop粒子，并且已按适应度排序
num=size(pop,2);  %pop为结构体
f_sum=0;
names=fieldnames(pop);
for i=1:num
    pop(i).x=getfield(pop(i),names{n});
end
for i=1:num
   f_sum=f_sum+pop(i).fitness;
end
f_worst=pop(n).fitness;
alpha=norm(UB-LB)/(f_worst/f_sum);
for i=1:num
    for j=1:num
        if i~=j
            reer(i,j)=alpha*(pop(i).fitness/f_sum)/norm(pop(j).x-pop(i).x);
            d(i,j)=norm(pop(j).x-pop(i).x);
        else
            reer(i,j)=0;
            d(i,j)=0;
        end
    end
end
k=0;
ser=1:num;
y1=[];
while 1
    if isempty(ser)
        break
    elseif length(ser)==1
        y1=[y1 ser];
        break
    end
    a=max(max(reer));     %最大的reer值为最优
    [b,a1]=find(reer==a);
    d_avg=sum(d(b,:))/(length(d(b,:))-1);
    temp1=find(d(b,:)<d_avg);
    if length(temp1)>2   %小生境群体
        k=k+1;
        y{k}=ser(temp1);
        reer=redu(reer,temp1,'rc');
        d=redu(d,temp1,'rc');
        ser=redu(ser,temp1,'c');
    else                   %一般群体
        y1=[y1 ser(temp1)];
        reer=redu(reer,temp1,'rc');
        d=redu(d,temp1,'rc');
        ser=redu(ser,temp1,'c');     
    end
end
for i=1:k
    for j=1:length(y{i})
        new{i}(j)=pop(y{i}(j));
    end
end
new1=pop(y1);

  


