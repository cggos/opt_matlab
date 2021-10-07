function [best_x,fval]=binfish(fun,fishnum,LB,UB,nvar,max_iterm,visual,step,try_number,delta)   %人工鱼群算法求TSP,n为维数
esp1=1e-8;
if isempty(visual)
    parameter.visual=ceil(unifrnd(nvar*0.4,nvar*0.6,1));
else
    parameter.visual=ceil(visual);
end
if isempty(step)
    parameter.step=3;
else
    parameter.step=ceil(step);
end
parameter.try_number=try_number;
parameter.delta=delta;
CodeLen=nvar*max(ceil(log2((UB-LB)./esp1 + 1)));   %二进制长度
for i=1:fishnum
   afish(i,:)=rand(1,CodeLen)<0.5;
   y(i)=fun(afish(i,:));
end
[best_val best_num]=max(y);
best_x=afish(best_num,:);      %最优鱼的路径
fval=best_val;
num=0;
for j=1:max_iterm
    for i=1:fishnum
       afish(i,:)=fishevaluate3_1(fun,afish(i,:),afish,parameter);
       y(i)=fun(afish(i,:));
    end
    [f,idex]=sort(y,'descend');
   if f(1)>fval
       best_x=afish(idex(1),:);
       fval=f(1);
       num=0;
   elseif abs(f(1)-fval)<=1e-4
       num=num+1;
   end
   if num==2
      lamda=rand(1,CodeLen);
      afish(idex(end),:)=(lamda.*afish(idex(end),:)+(1-lamda).*afish(idex(1),:))<0.5;
      y(idex(end))=fun(afish(idex(end),:));
      if y(idex(end))>fval
          fval=y(idex(end));
          best_x=afish(idex(end),:);
      end
      num=0;
   end
end
fval=fun(best_x);
best_x=Dec(LB,UB,best_x,CodeLen);


