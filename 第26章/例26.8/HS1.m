function [best_x,fval]=HS1(fun,popsize,iter_max,alpha,LB,UB)   %»ìºÏÍÜÌø¡¢Ï¸¾úÃÙÊ³ºÍÉùËÑË÷»ìºÏËã·¨
NC=size(LB,1);
for i=1:popsize     %¼ÇÒä¿â
    harmony(i,:)=LB'+(UB-LB)'.*rand(1,NC);
    fitness(i,1)=fun(harmony(i,:));
end
[a,b]=sort(fitness,'ascend');
best_x=harmony(b(1),:);
fval=a(1);
worst=harmony(b(end),:);
worstF=a(end);
harmony=redu(harmony,b(end),'r');
fitness=redu(fitness,b(end),'r');
for iter=1:iter_max(1)
    deltaG=1-(alpha(1)/alpha(2))^(1/iter);
    alphaG=(1-deltaG)*alpha(2);
    for j=1:iter_max(2)
       new=worst+alphaG*rand*(best_x-worst);
       new=boundtest(new,LB,UB);
       newF=fun(new);
       if newF<worstF
           if newF<fval
               best_x=new;
               fval=newF;
               worst=best_x;
               worstF=fval;
           else
               worst=new;
               worstF=newF;
           end
       else
          new=worst-alphaG*rand*(best_x-worst);
          new=boundtest(new,LB,UB);
          newF=fun(new); 
          if newF>worst
              worst=LB'+(UB-LB)'.*rand(1,NC);
              worstF=fun(worst);
          else
              if newF<fval
                 best_x=new;
                 fval=newF;
                 worst=best_x;
                 worstF=fval;
              else
                 worst=new;
                 worstF=newF;
              end
          end
       end
    end
    pop=[harmony;worst];
    y=[fitness;worstF];
    [a,b]=sort(y,'ascend');
    worst=pop(b(end),:);
    worstF=a(end);
    harmony=redu(pop,b(end),'r');
    fitness=redu(y,b(end),'r');
end
    
    
            
            
