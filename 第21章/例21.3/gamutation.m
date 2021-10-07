function pop=gamutation(pop,pm,varargin)
popsize=size(pop,2);
genelen=length(pop(1).x);
if strcmp(varargin{1},'b')
    for i=1:popsize
      for j=1:genelen
        test=rand;
        if test<pm
            pop(i).x(j)=abs(pop(i).x(j)-1);
        end
      end
    end
else
    LB=varargin{1};UB=varargin{2};
    iterm=varargin{3};iterm_max=varargin{4};
    for i=1:popsize
      for j=1:genelen
         if rand<pm 
            v1=pop(i).x(j)-UB(j);
            v2=LB(j)-pop(i).x(j);
            fg=rand*(1-iterm/iterm_max)^2;
            if rand>0.5
               pop(i).x(j)=pop(i).x(j)+v1*fg;
            else
               pop(i).x(j)=pop(i).x(j)+v2*fg;
            end
            pop(i).x(j)=boundtest(pop(i).x(j),LB(j),UB(j));
         end
      end
   end
end          