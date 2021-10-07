function [best_x,fval]=newfish(fun,fishnum,max_iterm,LB,UB,visual,try_number,delta)    %改进的鱼群算法
%prompt={'鱼群数';'最大迭代数';'LB';'UB';'视野';'步长';'试验次数';'拥挤度';};
%name='输入算法各参数';
%defaultanswer={'30','200','[]','[]','[]','[]','[]','[]'};
%answer=inputdlg(prompt,name,1,defaultanswer);
%fishnum=str2num(answer{1});
%max_iterm=str2num(answer{2});
%LB=str2num(answer{3});
%UB=str2num(answer{4});
%parameter.visual=str2num(answer{5});
%parameter.step=str2num(answer{6});
%parameter.try_number=str2num(answer{7});
%parameter.delta=str2num(answer{8});
parameter.visual=visual;
%parameter.step=step;
parameter.try_number=try_number;
parameter.delta=delta;
n=size(LB,1);    %维数
for i=1:fishnum
    afish(i,:)=LB'+(UB-LB)'.*rand(1,n);   %初始值
    y(i)=fun(afish(i,:)); 
end
[best_value,best_num]=max(y);
best_x=afish(best_num,:);
fval=best_value;         %全局极值
num=0;
k=1;
for j=1:max_iterm
  a=exp(-30*(j/max_iterm)^5);
  visual1=parameter.visual*a+0.001;
  % parameter.step=parameter.step*a+0.002;
   parameter.delta=0.95*parameter.delta;
  %parameter.step=(0.5+1/j)*parameter.step;
   for i=1:fishnum
       if ~isequal(afish(i,:),best_x)
           parameter.visual=fishdstc(afish(i,:),best_x);
       else
          parameter.visual=visual1; 
       end
       parameter.step=0.90*parameter.visual;
       afish(i,:)=fishevaluate(fun,afish(i,:),afish,LB,UB,best_x,parameter);
       y(i)=fun(afish(i,:));
       if y(i)>fval
           best_x=afish(i,:);
           fval=y(i);
      end
   end
   [f,idex]=sort(y,'descend');
   if f(1)>fval
       best_x=afish(idex(1),:);
       fval=f(1);
       num=0;
   else
       num=num+1;
   end
   if num==2
      k=0.99*k;
      afish(idex(end),:)=(1+k*rand).*afish(idex(end),:);    %高斯变异
      afish(idex(end),:)=boundtest(afish(idex(end),:),LB,UB);
      y(idex(end))=fun(afish(idex(end),:));
      if y(idex(end))>fval
          fval=y(idex(end));
          best_x=afish(idex(end),:);
      end
      num=0;
   end
   if j==max_iterm
       a=find(y<0.1*fval);
       for i=1:length(a)
           afish(a(i),:)=LB'+(UB-LB)'.*rand(1,n);   %重新初始化
           y(a(i))=fun(afish(a(i),:));
           if y(a(i))>fval
               fval=y(a(i));
               best_x=afish(a(i),:);
           end
       end
       %afish=redu(afish,a,'r');　　　%吞食
       %fishnum=fishnum-length(a);
       %y=redu(y,a,'c');
   end
end



