function [best_x,fval]=HS(fun,popsize,iter_max,LB,UB)   %∫Õ…˘À—À˜À„∑®
NC=size(LB,1);
HMCRmin=0.1;HMCRmax=0.5;
PARmin=0.9;PARmax=0.99;
BWmin=0.01;BWmax=0.1;
for i=1:popsize     %º«“‰ø‚
    harmony(i,:)=LB'+(UB-LB)'.*rand(1,NC);
    fitness(i,1)=fun(harmony(i,:));
end
[a,b]=sort(fitness,'ascend');
best_x=harmony(b(1),:);
fval=a(1);
for iter=1:iter_max
    PAR=PARmin+(PARmax-PARmin)*iter/iter_max;
    HMCR=HMCRmin+(HMCRmax-HMCRmin)*iter/iter_max;
    BW=BWmax*exp(iter*log(BWmin/BWmax)/iter_max);
    for i=1:popsize
        if rand<HMCR
            num=ceil(popsize*rand(1,NC));
            for j=1:NC
                new(i,j)=harmony(num(j),j);
            end
            if rand<PAR
                new(i,:)=new(i,:)+BW.*rand(1,NC);
            end
        else
            new(i,:)=LB'+(UB-LB)'.*rand(1,NC);
        end
        y1(i,1)=fun(new(i,:)); 
    end
    pop=[harmony;new];
    y=[fitness;y1];
    [a1,b1]=sort(y,'ascend');
    if a1(1)<fval
        best_x=pop(b1(1),:);
        fval=a1(1);
    end
    pop=pop(b1,:);
    y=y(b1,1);
    n=size(pop,1);
    harmony=redu(pop,popsize+1:n,'r');
    fitness=redu(y,popsize+1:n,'r');
end
    
    
            
            
