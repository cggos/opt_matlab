
function [best_x,fval]=BA(fun,popsize,iterm_max,LB,UB)     %òùòðËã·¨
NC=size(LB,1);
A=unifrnd(1,2);
alpa=unifrnd(0.8,0.95); 
gamma=unifrnd(0.8,0.95);
r=unifrnd(0,1);
fmin=0;fmax=2; 
for i=1:popsize
    bat(i).x=LB'+(UB-LB)'.*rand(1,NC);
    bat(i).fitness=fun(bat(i).x);
    bat(i).v=rand(1,NC);
    y(i)=bat(i).fitness;
end
[a,b]=min(y);
cbest=bat(b);
favg=mean(y);
for iterm=1:iterm_max
     for i=1:popsize
         if bat(i).fitness>favg
             omiga=0.9;
         else
             omiga=0.4-0.5*(bat(i).fitness-a)/(favg-a);
         end
         bat(i).f=fmin+(fmax-fmin)*rand;   
         bat(i).v=bat(i).v+(bat(i).x-cbest.x)*bat(i).f; 
         bat(i).x=bat(i).x+rand*bat(i).v; 
         bat(i).x=boundtest(bat(i).x,LB,UB);
         bat(i).fitness=fun(bat(i).x);
         y(i)=bat(i).fitness;
         if rand>r 
             new=cbest.x+A*unifrnd(-1,1,1,NC);
             new=boundtest(new,LB,UB);
             Fnew=fun(new);
             if Fnew<=bat(i).fitness&&rand<A 
                bat(i).x=new; 
                bat(i).fitness=Fnew;
                y(i)=Fnew;
                A=alpa.*A; 
                r=r.*(1-exp(-gamma.*iterm)); 
             end    
         end
     end
     [a,b]=min(y);
     if a<cbest.fitness
         cbest=bat(b);
     end
     favg=mean(y);
end
best_x=cbest.x;
fval=cbest.fitness;
  

