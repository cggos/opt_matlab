function xmin=newchoas(fun,x0,LB,UB,N1,N2,esp)     %牛顿－混沌法求解方程组 
NC=length(x0);
for i=1:N1
    for j=1:N2
      [fx dfx]=fun(x0);
      y=(dfx)\-fx;
      if norm(y)<esp
         xmin(i,:)=x0;
         break
      end
      x0=x0+y';
    end
    a=exist('xmin','var');
    if a==0
        xmin(i,:)=inf.*ones(1,NC);
        x0=LB'+(UB-LB)'.*rand(1,NC);
    end
    z=(x0-LB')./(UB-LB)';
    z=4*z.*(1-z);
    x0=LB'+z.*(UB-LB)';
end
xmin=delsample(xmin); 
    
    