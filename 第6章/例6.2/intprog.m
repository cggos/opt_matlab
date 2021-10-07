function [x,val,stats]=intprog(f,A,b,I,Aeq,Beq,lb,ub)  %整数规划或混合整数规划
%I为整数约束（即整数约束的变量序号）
if nargin<8
    ub=[];
    if nargin<7
       lb=[];
       if nargin<6
           Beq=[];
           if nargin<5
               Aeq=[];
               if nargin<4
                   I=1:length(f);
               end
           end
       end
    end
end
options=optimset('display','off');
[x0,val0,exitflag]=linprog(f,A,b,Aeq,Beq,lb,ub,[],options);
if exitflag<0
    disp('没有合适整数解');
    x=x0;
    val=val0;
    stats=exitflag;
    return;
else
    bound=inf;
    [x,val,stats]=branch(f,A,b,I,x0,val0,bound,Aeq,Beq,lb,ub);
end

function [newx,newval,stats,newbound]=branch(f,A,b,I,x,val,bound,Aeq,Beq,lb,ub)
options=optimset('display','off');
[x0,val0,stats0]=linprog(f,A,b,Aeq,Beq,lb,ub,[],options);
if stats0<0||val0>=bound
    newx=x;
    newval=val;
    newbound=bound;
    stats=stats0;
else
    index=myinteger(x0(I));
    if isempty(index)
       newx=x0;
       newval=val0;
       newbound=val0;
       stats=1;
    else
       n=I(index(1));
       add_A=zeros(1,length(f));
       add_A(n)=1;
       A=[A;add_A];
       b=[b;floor(x(n))];
       [x1,val1,stats1,bound1]=branch(f,A,b,I,x0,val0,bound,Aeq,Beq,lb,ub);
       A(end,:)=[];
       b(end,:)=[];
       stats=stats1;
       if stats1>0&&bound1<bound
           newx=x1;
           newval=val1;
           bound=val1;
           newbound=bound1;
       else
           newx=x0;
           newval=val0;
           newbound=bound;
       end
       A=[A;-add_A];
       b=[b;-ceil(x(n))];
       [x2,val2,stats2,bound2]=branch(f,A,b,I,x0,val0,bound,Aeq,Beq,lb,ub);
       if stats2>0&&bound2<bound
           stats=stats2;
           newx=x2;
           newval=val2;
           newbound=bound2;
       end
    end
end


