function [xmin,minf]=zoutendijk1(fun,A,b,gfun,x0,x_syms,esp1,esp2,type)   
%zoutendijk可行方向法,type=1线性约束,此时gfun=[];
%type=2非线性约束,此时需输入矩阵形式的约束函数的符号变量表达式gfun
syms x
if nargin==6
    if isempty(gfun)
       type=1;
       esp1=1e-6;
       esp2=1e-6;
    elseif ~isempty(gfun)
       type=2;
       esp1=1e-6;
       esp2=1e-6;
    end
elseif nargin==7
    if isempty(gfun)
       type=1;
       esp2=1e-6;
    elseif ~isempty(gfun)
       type=2;
       esp2=1e-6;
    end
elseif nargin==8
    if isempty(gfun)
       type=1;
    elseif ~isempty(gfun)
       type=2;
    end
end
z=myjacobian1(fun,x_syms);
n=length(x0);
k=0;
while k<500
    if type==1
       if size(A,2)~=size(x0,1)
           gi=A*x0'-b;
       else
           gi=A*x0-b; 
       end
       y1=findzeros(gi,0);
       g0=eval(subs(z,x_syms,x0))';
       if isempty(y1)
          if norm(g0)<esp1
             xmin=x0;
             minf=eval(subs(fun,x_syms,xmin));
             break
          else
             dk=-g0;
          end
       else
          A1=A(y1,:);
          b1=zeros(length(y1),1);
          A2=redu(A,y1,'r');   %不等于0
          b2=redu(b,y1,'r');
          L=-ones(n,1);u=ones(n,1);
          options=optimset('Display','off');
          dk=linprog(g0,-1*A1,b1,[],[],L,u,[],options);
          zk=eval(subs(z,x_syms,x0))';
          zk=zk'*dk;
          if abs(zk)<esp2
             xmin=x0;
             minf=eval(subs(fun,x_syms,xmin));
             break
          end
       end
       b_av=b2-A2*x0';
       d_av=A2*dk;
       id=find(d_av<0);
       d_av=d_av(id,:);
       b_av=b_av(id,:);
       alpha_av=min(b_av./d_av);
       alpha=mygoldcut1(fun,x0,dk,0,alpha_av,x_syms);
       xk=x0+alpha*dk'; 
    elseif type==2
       gi=eval(subs(gfun,x_syms,x0));
       g0=eval(subs(z,x_syms,x0))';
       y1=findzeros(gi,0);
       if isempty(y1)
          if norm(g0)<esp1
             xmin=x0;
             minf=eval(subs(fun,x_syms,xmin));
             break
          else
             dk=-g0;
          end
       else
          Igi=gfun(y1);
          zI=myjacobian1(Igi,x_syms);
          gI=eval(subs(zI,x_syms,x0));
          m1=length(y1);
          c=[zeros(n,1);1];
          A2=[g0' -1];
          A3=[-1*gI -1*ones(m1,1)];
          A4=[ones(1,n) 0];
          A5=[-1*ones(1,n) 0];
          b1=zeros(m1+1,1);
          b2=[b1;1;-1];
          options=optimset('Display','off');
          dk=linprog(c,[A2;A3;A4;A5],b2,[],[],[],[],[],options);
          zk=dk(end);
          dk=dk(1:end-1);
          if abs(zk)<esp2
              xmin=x0;
              minf=eval(subs(fun,x_syms,xmin));
              break     
          end
       end
       phi=subs(fun,x_syms,x0+x*dk');
       [a,b]=myJT1(phi,0,0.2);
       al=mygoldcut1(fun,x0,dk,a,b,x_syms); 
       xk=x0+al*dk';
    end
    x0=xk;
    k=k+1;
end

















    