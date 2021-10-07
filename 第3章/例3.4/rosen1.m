function [xmin,minf]=rosen1(fun,A,b,E,x0,x_syms,esp)  %Rosen投影法
if nargin==6
    esp=1e-6;
end
z=myjacobian1(fun,x_syms);
n=length(x0);
I=eye(n);
if size(A,2)~=size(x0,1)
    gi=A*x0'-b;
else
    gi=A*x0-b; 
end
y1=findzeros(gi,0);
g0=eval(subs(z,x_syms,x0))';
flag=1;
k=0;
while k<500
    if norm(g0)<esp
       xmin=x0;
       minf=eval(subs(fun,x_syms,xmin));
       break;
    end
    if ~isempty(y1)
       if flag==1
          A1=A(y1,:);
          A2=redu(A,y1,'r');   %不等于0
          b2=redu(b,y1,'r');
       elseif flag==0
          A1=redu(A1,y1(end),'r');
       end
       m=size(A1,1);
       if ~isempty(E)
          M=[A1;E];
          num=size(E,1);
       else
          M=A1;
          num=0;
       end
    else
       if isempty(E)
          M=[];
          num=0;
       else
          M=E;
          num=size(E,1);
       end
    end
    if isempty(M)
       P=I;
    else
       P=I-M'*(M*M')^-1*M;
    end
    dk=-P*g0;
    if norm(dk)>1e-3
       b_av=b2-A2*x0';
       d_av=A2*dk;
       id=find(d_av<0);
       d_av=d_av(id,:);
       b_av=b_av(id,:);
       alpha_av=min(b_av./d_av);
       alpha=mygoldcut1(fun,x0,dk,0,alpha_av,x_syms); 
       xk=x0+alpha*dk';
       x0=xk;
       k=k+1;
       flag=1;
       if size(A,2)~=size(x0,1)
          gi=A*x0'-b;
       else
          gi=A*x0-b; 
       end  
       y1=findzeros(gi,0);
       g0=eval(subs(z,x_syms,x0))';
    else
      omiga=(M*M')^-1*M*g0;
      if num==0
          lamda=omiga;
      else
          lamda=redu(omiga,m+1:m+num,'r');
      end
      [y1,A4]=find(lamda<0);
      if isempty(y1)
         xmin=x0;
         minf=eval(subs(fun,x_syms,xmin));
         break;
      else
         flag=0;
      end
   end
end

    