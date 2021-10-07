function Candidate=selectCandidate(v,Candidate_Num,varnum,iterm,iterm_max,LB,UB,pm,type)
for i=1:Candidate_Num
    Candidate(i).x=v;
    if type==1
       for j=1:varnum
          if rand>pm
            v1=v(j)-UB(j);
            v2=LB(j)-v(j);
            fg=rand*(1-iterm/iterm_max)^2;
            if rand>0.5
               Candidate(i).x(j)=v(j)+v1*fg;
            else
               Candidate(i).x(j)=v(j)+v2*fg;
            end
          else
             Candidate(i).x(j)=v(j)+rand*(UB(j)-LB(j));
          end
          if Candidate(i).x(j)>=UB(j)
             Candidate(i).x(j)=LB(j)+(Candidate(i).x(j)-UB(j));
          elseif Candidate(i).x(j)<=LB(j)
             Candidate(i).x(j)=UB(j)-(LB(j)-Candidate(i).x(j));
          end
       end
    elseif type==2
        b=ceil(rand*varnum);
        Candidate(i).x(b)=LB(b)+rand*(UB(b)-LB(b));
        if Candidate(i).x(b)>=UB(b)
              Candidate(i).x(b)=LB(b)+(Candidate(i).x(b)-UB(b));
        elseif Candidate(i).x(b)<=LB(b)
             Candidate(i).x(b)=UB(b)-(LB(b)-Candidate(i).x(b));
        end
    end
end



