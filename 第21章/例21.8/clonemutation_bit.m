function cloneantibody=clonemutation_bit(cloneantibody,pm)
cloneantibodynum=size(cloneantibody,2);
NC=length(cloneantibody(1).x);
for i=1:cloneantibodynum
    for j=1:NC
        p=rand;
        if p<pm
           cloneantibody(i).x(j)=abs(cloneantibody(i).x(j)-1);
        end
    end
end
