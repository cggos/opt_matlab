function [bestx,fval]=SMSA(fun,x0,LB,UB)    %单纯形法－模拟退火算法
if isempty(x0)
    NC=size(LB,1);
else
    NC=length(x0(1,:));
end
if size(x0,1)~=NC+1&&isempty(LB)
    s=eye(NC);
    for i=1:NC+1
        if i==1
           x(i,:)=x0;
        else
           x(i,:)=x0+nuifnd(-1,1)*s(i-1,:);
           x(i,:)=boundtest(x(i,:),LB,UB);
        end
        fit(i,1)=fun(x(i,:));
    end
elseif isempty(LB)
    for i=1:NC+1
        x(i,:)=x0(i,:);
        fit(i,1)=fun(x(i,:));
    end
elseif isempty(x0)
    s=eye(NC);
    x(1,:)=LB'+(UB-LB)'.*rand(1,NC);
   for i=2:NC+1
       x(i,:)=x(1,:)+unifrnd(-1,1)*s(i-1,:);
       x(i,:)=boundtest(x(i,:),LB,UB);
       fit(i,1)=fun(x(i,:));
   end
end
t=1;
for iter=1:10000
    if t<0.1^30
        break;
    end
    for i=1:150
       [a,b]=sort(fit,'descend');
       fh=a(1)-t*log(rand);fl=a(end)-t*log(t);
       xh=x(b(1),:);x_retain=x(b(2:end),:);f_retain=a(2:end);
       xl=(x(b(end),:));
       xc=(sum(x)-xh)./NC;
       gama=1;lamda=0.75;mu=2;   %单纯形法参数
       xr=2*xc-gama*xh;
       xr=boundtest(xr,LB,UB);
       fr=fun(xr)+t*log(rand);
       if fr>=fh
           xs=(1-lamda)*xh+lamda*xr;
           fs=fun(xs)+t*log(rand);
       else
          xe=xh+mu*(xr-xh);
          xe=boundtest(xe,LB,UB);
          fe=fun(xe)+t*log(rand);
          if fe<=fr
             xs=xe;
             fs=fe;
          else
             xs=xr;
             fs=fr;
          end
       end
       if fs<fh
          xh=xs;
          fh=fun(xh);
          x=[xh;x_retain];
          fit=[fh;f_retain];
      else
         e=0;
         for j=1:NC+1
            e=e+(fit(j,1)-fl)^2;
         end
         if sqrt(e/(NC+1))<=1e-6
            break
         else
            for j=1:NC+1
                x(j,:)=(x(j,:)+xl)/2;
                fit(j,1)=fun(x(j,:));
            end
         end         
       end
    end
    for i=1:100
        w=unifrnd(0,1,1,NC);
        u=unifrnd(-1,1,1,NC);
        totalW=sqrt(sum(w.^2));
        for j=1:NC
            new(j)=xl(j)+w(j)*t*(1/u(j)^3-1)/totalW;
        end
        new=boundtest(new,LB,UB);
        newF=fun(new);
        df=newF-fl;
        if df<0
           xl=new;
           fl=newF;
        elseif min(1,exp(-df/t))>rand
           xl=new;
           fl=newF;
        end
    end
    t=0.95*t;
    x(1,:)=xl;
    for i=2:NC+1
        x(i,:)=x(1,:)+unifrnd(-1,1)*s(i-1,:);
        x(i,:)=boundtest(x(i,:),LB,UB);
        fit(i,1)=fun(x(i,:));
    end
end
bestx=xl;fval=fun(xl);
        

    



