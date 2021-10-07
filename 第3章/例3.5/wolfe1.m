function [xmin,minf]=wolfe1(fun,A,x0,x_syms,esp)
if nargin==4
    esp=1e-4;
end
z=myjacobian1(fun,x_syms);
n=length(x0);
m=size(A,1);
I=find(x0>0);
k=0;
while k<500
   I1=redu(1:n,I,'c');
   xB=x0(I);xN=redu(x0,I,'c');
   B=A(:,I);N=redu(A,I,'c');
   g0=eval(subs(z,x_syms,x0))';
   gB=g0(I,:);gN=redu(g0,I,'r');
   r=gN-(B^-1*N)'*gB;
   for i=1:n-m
      if r(i)>0
         dN(i,1)=-xN(i)*r(i);
      else
         dN(i,1)=-r(i);
      end
   end
   dB=-B^-1*N*dN;
   dk(I,1)=dB;dk(I1,1)=dN;
   if norm(dk)<esp
       xmin=x0;
       minf=eval(subs(fun,x_syms,xmin));
       break
   else
       a1=find(dB<0);a2=find(dB>0);
       if ~isempty(a1)
          b2=min(-xB(a1)'./dB(a1));
       else
           b2=inf;
       end
       if ~isempty(a2)
          b3=max(-xB(a2)'./dB(a2));
       else
           b3=inf;
       end
       L1=min(b2,b3);
       id1=find(dB<0);id2=find(dN<0);
       dkj=[];xkj=[];
       if ~isempty(id1)
          dkj=[dkj dB(id1)'];
          xkj=[xkj xB(id1)];
       end
       if ~isempty(id2)
          dkj=[dkj dN(id2)'];
          xkj=[xkj xN(id2)];
       end
       alpha_av=min(-xkj./dkj);
       alpha=mygoldcut1(fun,x0,dk,0,alpha_av,x_syms);
       xk=x0+alpha*dk';
       if abs(alpha-L1)<1e-6
          a3=findzeros(xk,0);
          y3=mycompare1(a3,I);
          y4=[];
          for i=1:length(y3)
             y4=[y4 find(I==y3(i))];
          end
          I=redu(I,y4,'c');
          [a4,b4]=max(xk(I1));
          I=sort([I b4]);
       end
       x0=xk;
       k=k+1;
   end
end




