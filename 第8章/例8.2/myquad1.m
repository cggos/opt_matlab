function [xk,lamda,minf]=myquad1(H,c,Ae,be,Ai,bi,x0)   %二次规划的有效集方法
ne=length(be);
ni=length(bi);
I=ones(ni,1);
for i=1:ni
    if Ai(i,:)*x0-bi(i)>1e-6   
        I(i)=0;
    end
end
k=0;
while k<500
    Aee=[];
    if ne>0
        Aee=Ae;
    end
    for j=1:ni
        if I(j)>0    %有效集为
            Aee=[Aee;Ai(j,:)];
        end
    end
    gk=H*x0+c;
    [m1,n1]=size(Aee);
    [dk,lamk]=subquad(H,gk,Aee,zeros(m1,1));
    if norm(dk)<1e-6
        y=0;
        if length(lamk)>ne
            [y,jk]=min(lamk(ne+1:length(lamk)));
        end
        if y>=0
            flag=0;
        else
            flag=1;
            for i=1:ni
                if I(i)&&ne+sum(I(1:i))==jk
                    I(i)=0;
                    break;
                end
            end
        end
        k=k+1;
    else
        flag=1;
        alpha=1.0;tm=1.0;
        for i=1:ni
            if I(i)==0&&Ai(i,:)*dk<0
                tm1=(bi(i)-Ai(i,:)*x0)/(Ai(i,:)*dk);
                if tm1<tm
                    tm=tm1;
                    ti=i;
                end
            end
        end
        alpha=min(alpha,tm);
        xk=x0+alpha*dk;
        x0=xk;
        if tm<1
            I(ti)=1;
        end
    end
    if flag==0
        break
    end
    k=k+1;
end
lamda=[lamk(1:ne);zeros(ni,1)];
p=find(I>0);
s=size(p,1);
for i=1:s
    lamda(ne+p(i))=lamk(ne+i);
end
minf=0.5*xk'*H*xk+c'*xk;


function [x,lambda]=subquad(H,c,Ae,be)
ginvH=pinv(H);
[m,n]=size(Ae);
if m>0
    rb=Ae*ginvH*c+be;
    lambda=pinv(Ae*ginvH*Ae')*rb;
    x=ginvH*(Ae'*lambda-c);
else
    x=-ginvH*c;
    lambda=0;
end

        
                
                    
            
        