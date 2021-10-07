function y=createCOA(L,n,type)     %产生n维的混沌序列
if strcmp(type,'logistic')
    y=rand(1,n);
    for i=1:L
      y=4.*y.*(1-y);
    end
elseif strcmp(type,'tent')
    y=unifrnd(0.5,1,1,n);
    for i=1:L
      y=2.*(1-y);
    end
elseif strcmp(type,'henon')
    x=zeros(1,n);y=zeros(1,n);a=unifrnd(0,1.4,1,n);
    for i=1:L
        xm=x;
        ym=y;
        x=ym+1-a.*xm.*xm;
        y=0.3.*xm;
    end
elseif strcmp(type,'zhb')    %账篷
    a=rand(1,n);y=0.34.*ones(1,n);
    for i=1:L
       y=a-(1+a).*abs(y);
    end
elseif strcmp(type,'kent')
    a=unifrnd(0,0.5,1,n);y=rand(1,n);
    for i=1:L
       for j=1:n
          if y(j)>0&&y(j)<=a(j)
              y(j)=y(j)/a(j);
          elseif y(j)>a(j)&&y(j)<1
              y(j)=(1-y(j))/(1-a(j));
          end
       end
    end
elseif strcmp(type,'sin')
    y=unifrnd(-1,1,1,n);
    for i=1:L
      y=sin(5.56./y);
    end
end
    
        