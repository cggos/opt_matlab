function minx=multgoal(f,A,b,I)      %多目标规划的单纯形法,f是目标函数,I为虚拟变量的序号
if nargin==3
    I=[];
end
n=size(f,1);    %求目标函数的层数
[r,c]=size(A);
minx=zeros(1,c);
m=findeye(A);
for i=1:n
  cI=f(i,m(:,1));
  for j=1:c
     sigma(i,j)=cI*A(:,j)-f(i,j);
  end
end
for i=1:n
    while 1
      [a,id]=find(sigma(i,:)>0);
      if isempty(a)
          break
      end
      b1=find(A(:,id(1))>0);
      [a4,ser]=redu(1:r,b1,'c');
      [a,a3]=sort(b(b1)./A(b1,id(1)));
      b3=find(a==a(1));
      idex=ser(a3(b3(end)));
      y=progturn(A,b,sigma(i,:),[],idex,id(1));
      A=y(1:end-1,1:c);
      b=y(1:end-1,c+1);
      m=findeye(A);
      for j=1:n
         cI=f(j,m(:,1));
         for k=1:c
            sigma(j,k)=cI*A(:,k)-f(j,k);
         end
      end
    end
end
y=findeye(A);
minx(y(:,1)')=b(y(:,2));
minx=redu(minx,I,'c');



    
    