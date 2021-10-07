function pop=gacrossover(pop,px,xtype,type,LB,UB)
if nargin==4
    LB=[];UB=[];
end
popsize=size(pop,2);
ndim=length(pop(1).x);
halfpop=popsize/2;
if strcmp(type,'b')
  if xtype==1   %单点交换
    randlist=rand((halfpop),1);
    for i=1:halfpop
        a=(i*2)-1;
        xpo=a+1;
        new(a).x(1:ndim)=pop(a).x(1:ndim);
        new(xpo).x(1:ndim)=pop(xpo).x(1:ndim);
        if (randlist(i)<px)
            xpoint=round((rand*ndim)+0.5);
            new(xpo).x(1:xpoint)=pop(a).x(1:xpoint);
            new(a).x(1:xpoint)=pop(xpo).x(1:xpoint);
        end
    end
  end
  if xtype==2     %两点交换
    randlist=rand((halfpop),1);
    for i=1:halfpop
        a=(i*2)-1;
        xpo=a+1;
        new(a).x(1:ndim)=pop(a).x(1:ndim);
        new(xpo).x(1:ndim)=pop(xpo).x(1:ndim);
        if (randlist(i)<px)
            xpoint=sort(round((rand(1,2)*ndim)+0.5));
            new(xpo).x(xpoint(1):xpoint(2))=pop(a).x(xpoint(1):xpoint(2));
            new(a).x(xpoint(1):xpoint(2))=pop(xpo).x(xpoint(1):xpoint(2));
        end
    end
  end
  if xtype==3
      for i=1:halfpop
        a=(i*2)-1;
        xpo=a+1;
        new(a).x(1:ndim)=pop(a).x(1:ndim);
        new(xpo).x(1:ndim)=pop(xpo).x(1:ndim);
        for j=1:ndim
            test=rand;
            if test<px
                new(xpo).x(j)=pop(a).x(j);
                new(a).x(j)=pop(xpo).x(j);
            end
        end
      end
  end
  pop=new;
elseif strcmp(type,'r')     %实数编码
    for i=1:popsize
        index=ceil(rand(1,2).*popsize);
        if rand<px
           pos=ceil(rand.*ndim);
           v1=pop(index(1)).x(pos);
           v2=pop(index(2)).x(pos);
           pop(index(1)).x(pos)=rand*v2+(1-rand)*v1;
           pop(index(2)).x(pos)=rand*v1+(1-rand)*v2;
           pop(index(1)).x(pos)=boundtest(pop(index(1)).x(pos),LB(pos),UB(pos));
           pop(index(2)).x(pos)=boundtest(pop(index(1)).x(pos),LB(pos),UB(pos)); 
        end  
    end
end

    
             
            
    