function cloneantibody=clonemutation(cloneantibody,iterm,iterm_max,pm,LB,UB)
cloneantibodynum=size(cloneantibody,2);
NC=length(cloneantibody(1).x);
for i=1:cloneantibodynum
    for j=1:NC
        p=rand;
        if p<pm
           v1=cloneantibody(i).x(j)-UB(j);
           v2=LB(j)-cloneantibody(i).x(j);
           fg=rand*(1-iterm/iterm_max)^2;
           if rand>0.5
              cloneantibody(i).x(j)=cloneantibody(i).x(j)+v1*fg;
           else
              cloneantibody(i).x(j)=cloneantibody(i).x(j)+v2*fg;
           end
           cloneantibody(i).x(j)=boundtest(cloneantibody(i).x(j),LB(j),UB(j));
        end
    end
end
