function pop=chaos_serch(fun,pop,LB,UB)  %»ìãçËÑË÷
popsize=size(pop,2);
for k=1:popsize
    fval=pop(k).fitness;
    best_x=pop(k).x;
    z=(best_x-LB')./(UB-LB)';%±ä»»µ½[0,1]Çø¼ä
    for i=1:1000  %´ÖËÑË÷
       z=4*z.*(1-z);  %²úÉún1Î¬µÄ»ìãç±äÁ¿
       mx=LB'+(UB-LB)'.*z;
       y=fun(mx);
       if y>fval
          best_x=mx;
          fval=y;
       end 
    end
    for i=1:300     %Ï¸ËÑË÷
        alpha=1-((i-1)/i)^2;
        x=best_x+alpha.*(z-0.5);
        y=fun(x);
        if y>fval
           best_x=x;
           fval=y;
        end
        z=4*z.*(1-z);
    end
    pop(k).x=best_x;
    pop(k).fitness=y;      
end
