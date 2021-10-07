function afish=fishswarm(fun,afish0,afish1,LB,UB,best_x,parameter)     %聚群,afish1为整个鱼群
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
      afish=afish0+r_step.*rand(1,m).*(afish_center-afish0+best_x-afish0)/norm(afish_center-2.*afish0+best_x);
      afish=boundtest(afish,LB,UB);
      return  
    end
end
afish=fishprey(fun,afish0,LB,UB,best_x,parameter);