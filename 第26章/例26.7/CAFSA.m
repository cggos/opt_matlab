function [best_x,fval]=CAFSA(fun,fishnum,max_iterm,visual,step,try_number,delta,LB,UB)
parameter.visual=visual;
parameter.step=step;
parameter.try_number=try_number;
parameter.delta=delta;
NC=size(LB,1);    %维数
z=rand(1,NC);
for i=1:fishnum
    z=4*z.*(1-z);
    afish(i,:)=LB'+(UB-LB)'.*z;
    y(i)=fun(afish(i,:));
end
[best_value,best_num]=max(y);
best_x=afish(best_num,:);
fval=best_value;         %全局极值
gama=unifrnd(0,0.5);
alpha=0.02;
Pb=0.95;
for j=1:max_iterm
    for i=1:fishnum
       afish(i,:)=fishevaluate6(fun,afish(i,:),afish,LB,UB,parameter);
       y(i)=fun(afish(i,:));
       if y(i)>fval
           best_x=afish(i,:);
           fval=y(i);
       elseif rand<Pb
           afish(i,:)=fishmove(afish(i,:),parameter,LB,UB);
           y(i)=fun(afish(i,:));
           if y(i)>fval
              best_x=afish(i,:);
              fval=y(i);
           end
       else
           afish(i,:)=afish(i,:)+rand*parameter.step*(best_x-afish(i,:))/norm(best_x-afish(i,:));
           afish(i,:)=boundtest(afish(i,:),LB,UB);
           y(i)=fun(afish(i,:));
           if y(i)>fval
              best_x=afish(i,:);
              fval=y(i);
           end
       end
    end
    Pb=0.98*Pb; 
    for i=1:30 
        LB1=LB;UB1=UB;
        for k=1:NC
            LBi(k,1)=best_x(k)-gama*(UB1(k,1)-LB1(k,1));
            if LBi(k,1)<LB1(k,1)
               LBi(k,1)=LB1(k,1);
            end
            UBi(k,1)=best_x(k)+gama*(UB1(k,1)-LB1(k,1));
            if UBi(k,1)<UB1(k,1)
               UBi(k,1)=UB1(k,1);
            end
        end
        LB1=LBi;UB1=UBi;
        x=(best_x-LB1')./(UB1-LB1)';
        for k=1:500
            y1=(1-alpha).*x+alpha.*z;
            z1=LB1'+(UB1-LB1)'.*y1;           
            y2=fun(z1);
           if y2<fval
              best_x=z1;
              x=(best_x-LB1')./(UB1-LB1)';
              fval=y2;
           end
           z=4*z.*(1-z);
        end
        alpha=alpha/2;
    end
end




function afish=fishevaluate6(fun,afish0,afish1,LB,UB,parameter)   %鱼的评价
afish1=fishfollow6(fun,afish0,afish1,LB,UB,parameter);
afish2=fishswarm6(fun,afish0,afish1,LB,UB,parameter);
afish3=fishprey6(fun,afish0,LB,UB,parameter);
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
    afish=afish0;
end

function afish=fishswarm6(fun,afish0,afish1,LB,UB,parameter)     %聚群,afish1为整个鱼群
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
     %r_step=abs(1-(fun(afish0)/fun(afish_center)))*parameter.step;
      r_step=parameter.step;
      afish=afish0+r_step.*rand(1,m).*(afish_center-afish0)/norm(afish_center-afish0);
      afish=boundtest(afish,LB,UB);
      return  
    end
end
afish=fishprey6(fun,afish0,LB,UB,parameter);

function afish=fishprey6(fun,afish0,LB,UB,parameter)    %觅食
n=length(afish0);
for i=1:parameter.try_number
    afish_next=afish0+parameter.visual.*rand(1,n);
    yj=fun(afish_next);
    yi=fun(afish0);
    if yi<yj
      %r_step=abs(1-yi/yj)*parameter.step;
       r_step=parameter.step;
       afish=afish0+r_step.*rand(1,n).*(afish_next-afish0)/norm(afish_next-afish0);
       afish=boundtest(afish,LB,UB);
       return
    end
end
afish=afish0;

function afish=fishfollow6(fun,afish0,afish1,LB,UB,parameter)   %追尾函数
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
    %r_step=abs(1-(fun(afish0)/f_max))*parameter.step;
    r_step=parameter.step;
    afish=afish0+r_step.*rand(1,m).*(afish1(max_i,:)-afish0)/norm(afish1(max_i,:)-afish0);
    afish=boundtest(afish,LB,UB);
    return
end
afish=afish0;


    
    
   