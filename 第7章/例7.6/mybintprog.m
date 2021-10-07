function [xmin,val]=mybintprog(f,A,b,x0)   %0-1规划的隐枚举法,
%增加过滤条件，即求函数的上界
[r,c]=size(A);
flag=0;
if nargin==3
   [a1,b1]=sort(f,'ascend');
   x0=zeros(c,1);
   for i=1:r
      x0(b1(i),1)=1;
      if A*x0<=b
         flag=1;
         break;
      else
         x0(b1(i),1)=0;
      end
   end
end
if nargin==4
   flag=1;
   if size(x0,1)==1
      x0=x0';
   end
end
if flag==0    %考虑找不到初始点
   minf=inf;
else
   minf=f*x0;
   xmin=x0;
end
val=minf;
if ~isempty(x0)
   A=[A;f];
   b=[b;minf];
end
for i=0:2^c-1
    x1=mydec2bin(i,c);
    if A*x1<=b
        f2=f*x1;
        if f2<minf
            minf=f2;
            if ~~isempty(x0)
               b(r+1,1)=minf;
            end
            xmin=x1';
            val=minf;
        else
            continue
        end
    else
        continue
    end
end

function y=mydec2bin(x,n)
s=dec2bin(x,n);
for i=1:n
    y(i)=str2double(s(i));
end
y=transpose(y);
  
        
        
        
    
  
    


