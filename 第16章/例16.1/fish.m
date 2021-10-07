function [best_x,fval]=fish(fun)    %鱼群算法
prompt={'鱼群数';'最大迭代数';'LB';'UB';'视野';'步长';'试验次数';'拥挤度';};
name='输入算法各参数';
defaultanswer={'30','200','[]','[]','[]','[]','[]','[]'};
answer=inputdlg(prompt,name,1,defaultanswer);
fishnum=str2num(answer{1});
max_iterm=str2num(answer{2});
LB=str2num(answer{3});
UB=str2num(answer{4});
parameter.visual=str2num(answer{5});
parameter.step=str2num(answer{6});
parameter.try_number=str2num(answer{7});
parameter.delta=str2num(answer{8});
n=size(LB,1);    %维数
for i=1:fishnum
    afish(i,:)=LB'+(UB-LB)'.*rand(1,n);   %初始值
    y(i)=fun(afish(i,:)); 
end
[best_value,best_num]=max(y);
best_x=afish(best_num,:);
fval=best_value;         %全局极值
for j=1:max_iterm
    for i=1:fishnum
       afish(i,:)=fishevaluate(fun,afish(i,:),afish,LB,UB,parameter);
       y(i)=fun(afish(i,:));
       if y(i)>fval
           best_x=afish(i,:);
           fval=y(i);
       end
   end
end

function afish=fishevaluate(fun,afish0,afish1,LB,UB,parameter)
afish1=fishfollow(fun,afish0,afish1,LB,UB,parameter);
afish2=fishswarm(fun,afish0,afish1,LB,UB,parameter);
afish3=fishprey(fun,afish0,LB,UB,parameter);
af_best=afish1;
if fun(afish2)>fun(af_best)
    af_best=afish2;
end
if fun(afish3)>fun(af_best)
    af_best=afish3;
end
if fun(af_best)>fun(afish0)
    afish=af_best;
else
   afish=fishmove(afish0,parameter,LB,UB);
end

function afish=fishprey(fun,afish0,LB,UB,parameter)    %觅食
n=length(afish0);
for i=1:parameter.try_number
    afish_next=afish0+parameter.visual.*rand(1,n);
    yj=fun(afish_next);
    yi=fun(afish0);
    if yi<yj
       r_step=abs(1-yi/yj)*parameter.step;
       afish=afish0+r_step.*rand(1,n).*(afish_next-afish0)/norm(afish_next-afish0);
       afish=boundtest(afish,LB,UB);
       return
    end
end
afish=fishmove(afish0,parameter,LB,UB);

function afish=fishswarm(fun,afish0,afish1,LB,UB,parameter)     %聚群,afish1为整个鱼群
[afishnum,m]=size(afish1);
n=0;
afish_center=0;
for i=1:afishnum
    if ~isequal(afish0,afish1(i,:))
       if fishdstc(afish1(i,:),afish0)<parameter.visual
          n=n+1;
          afish_center=afish_center+afish1(i,:);
       end
    end
end
if n~=0
   afish_center=afish_center/n;
   if (fun(afish_center)>fun(afish0)*parameter.delta*n)
      r_step=abs(1-(fun(afish0)/fun(afish_center)))*parameter.step;
      afish=afish0+r_step.*rand(1,m).*(afish_center-afish0)/norm(afish_center-afish0);
      afish=boundtest(afish,LB,UB);
      return  
    end
end
afish=fishprey(fun,afish0,LB,UB,parameter);

function afish=fishfollow(fun,afish0,afish1,LB,UB,parameter)   %追尾函数
[fishnum,m]=size(afish1);
n=0;
f_max=-inf;
max_i=1;
for i=1:fishnum
    if ~isequal(afish1(i,:),afish0)
      if (fishdstc(afish1(i,:),afish0)<parameter.visual)
        n=n+1;
        if fun(afish1(i,:))>f_max
            f_max=fun(afish1(i,:));
            max_i=i;
        end
      end
    end
end
if ((f_max/n)>(parameter.delta*fun(afish0)))
    r_step=abs(1-(fun(afish0)/f_max))*parameter.step;
    afish=afish0+r_step.*rand(1,m).*(afish1(max_i,:)-afish0)/norm(afish1(max_i,:)-afish0);
    afish=boundtest(afish,LB,UB);
    return
end
afish=fishprey(fun,afish0,LB,UB,parameter);

function afish=fishmove(afish0,parameter,LB,UB)      %随机
n=length(afish0);
afish=afish0+parameter.visual.*unifrnd(-1,1,1,n);
afish=boundtest(afish,LB,UB);

function out=fishdstc(af1,af2)
out=norm(af1-af2);












