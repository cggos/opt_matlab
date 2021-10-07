function [x,val]=gomory(A,b,c)     %割平面法求整数规划
[minx,val,A2,b1]=mysimplex(A,b,c);
while 1
   n=length(b1);
   y=zeros(1,n);
   for i=1:n
     if abs(round(b1(i))-b1(i))>1.0e-7
       y(i)=1;
     end
   end
   a=find(y==1);
   if ~isempty(a)     %不是整数解
      A=A2(1:end-1,1:end-1);
      [r,c1]=size(A);
      b=A2(1:end-1,end);
      m=findeye(A);
      [fk,id]=max(mydec(b(a)));
      id=a(id);
      nbase=redu(1:c1,m(:,2),'c');    %非基变量
      B=A(:,nbase);
      nba_num=length(nbase);
      y1=zeros(1,nba_num);
      for i=1:nba_num
         y1(i)=-mydec(B(id,i));
      end
      A=[A zeros(r,1)];
      A1=zeros(1,r+1);
      A1(nbase)=y1;
      A=[A;A1 1];
      b=[b;-fk];
      c=[c 0];
      [minx,b1]=dualsimplex(A,b,c);
      A2=minx.A;
      val=minx.f;
   else
      break
   end
end
x=minx;
      
   
   
   