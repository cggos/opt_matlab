function [best_x,fval]=FA(fun,varargin)   %Ó©»ð³æËã·¨
type=varargin{end};
popsize=varargin{1};
iterm_max=varargin{2};
LB=varargin{3};
UB=varargin{4};
if type==1
   rs=varargin{5};
   nt=5;gama1=0.6;beta=0.08;ro=0.4;step=0.3;
   rd=rs.*ones(popsize,1);
elseif type==2
    beta0=0.2;gama=1; 
end
NC=size(LB,1);l0=5;
for i=1:popsize
    fireworm(i).x=LB'+(UB-LB)'.*rand(1,NC);
    fireworm(i).fitness=fun(fireworm(i).x);
    if type==1
        fireworm(i).l=(1-ro)*l0+gama1*fireworm(i).fitness;
    elseif type==2
        fireworm(i).I=fireworm(i).fitness;
    end
    y(i)=fireworm(i).fitness;
end
[a,b]=sort(y,'descend');
cbest=fireworm(b(1));
favg=mean(y);
fireworm=fireworm(b);
if type==1
    for iterm=1:iterm_max
        step=step*(1-iterm/iterm_max)+10^(-4);
        for i=1:popsize
           chome=[];k=1;
           for j=1:popsize
               if j~=i
                  if (norm(fireworm(j).x-fireworm(i).x)<rd(i))&&(fireworm(j).l>fireworm(i).l) 
                      chome=[chome j];
                      p(k)=fireworm(j).l-fireworm(i).l;
                      k=k+1;
                  end
               end
           end
           if isempty(chome)
               fireworm(i).x=fireworm(i).x+unifrnd(-1,1,1,NC);
           else
               p=p./sum(p);
               for j=1:length(chome)   
                   pick=rand;
                   while pick==0    
                      pick=rand;        
                   end
                   pick=pick-p(j);        
                  if pick<0        
                      index=chome(j);
                      break;  
                  end 
               end
               if exist('index','var')==0
                  [a1,b2]=max(p);
                  index=chome(b2);
               end
               fireworm(index).x=localsearch(fun,fireworm(index).x,LB,UB,1);  %¾Ö²¿ËÑË÷
               fireworm(i).x=fireworm(i).x+step*(fireworm(index).x-fireworm(i).x)/(norm(fireworm(index).x-fireworm(i).x));
           end
           fireworm(i).x=boundtest(fireworm(i).x,LB,UB,1);
           fireworm(i).fitness=fun(fireworm(i).x);
           fireworm(i).l=(1-ro)*fireworm(i).l+gama1*fireworm(i).fitness;
           rd(i)=min(rs,max(0,rd(i)+beta*(nt-length(chome))));
           y(i)=fireworm(i).fitness;
        end
        [a,b]=sort(y,'descend');
        if a(1)>cbest.fitness
           cbest=fireworm(b(1));
        end
    end
elseif type==2
    for iterm=1:iterm_max
        alpha=0.25*0.95^iterm;
        omiga=exp(-abs(cbest.fitness/(cbest.fitness-favg)));
        for i=1:popsize
            if i==popsize
               beta=beta0*exp(-gama*norm(cbest.x-fireworm(i).x))*exprnd(1); 
               fireworm(i).x=omiga.*fireworm(i).x+beta.*(cbest.x-fireworm(i).x)+alpha.*rand(1,NC); 
            else
               x1=cbest.x+fireworm(i).x-fireworm(i+1).x;
               x=0.25.*(x1-fireworm(i).x)+fireworm(i).x;
               beta=beta0*exp(-gama*norm(x-fireworm(i).x))*exprnd(1);
               fireworm(i).x=omiga.*fireworm(i).x+beta.*(x-fireworm(i).x)+alpha.*rand(1,NC);
            end
            fireworm(i).x=boundtest(fireworm(i).x,LB,UB,1);
            fireworm(i).fitness=fun(fireworm(i).x);
            y(i)=fireworm(i).fitness;
        end
        [a,b]=sort(y,'descend');
        if a(1)>cbest.fitness
           cbest=fireworm(b(1));
        end
        favg=mean(y);
        fireworm=fireworm(b);
    end
end
best_x=cbest.x;
fval=cbest.fitness;
           
    
    
    
    