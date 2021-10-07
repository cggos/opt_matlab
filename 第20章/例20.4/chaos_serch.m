function pop=chaos_serch(fun,pop,LB,UB)  %��������
popsize=size(pop,2);
for k=1:popsize
    fval=pop(k).fitness;
    best_x=pop(k).x;
    z=(best_x-LB')./(UB-LB)';%�任��[0,1]����
    for i=1:1000  %������
       z=4*z.*(1-z);  %����n1ά�Ļ������
       mx=LB'+(UB-LB)'.*z;
       y=fun(mx);
       if y>fval
          best_x=mx;
          fval=y;
       end 
    end
    for i=1:300     %ϸ����
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
