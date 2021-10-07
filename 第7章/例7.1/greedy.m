function [y,val]=greedy(f,A,b,I)    %背包问题的贪心算法，
%%　f,A分别为价值和体积，b为约束条件,I是已确定的变量序号及对应的取值
if nargin==3
    I=[];
    n=length(f);
    ratio=f./A;
    a1=[];
end
if nargin==4
    if ~isempty(I)
       a1=find(I(:,2)==0);   %找取值为零的变量
       a2=find(I(:,2)==1);   %取值为1的变量
       n=length(f);
       for i=1:n
         if i==mycompare1(i,I(a1,1))
            ratio(i)=0;
         elseif i==mycompare1(i,I(a2,1))
            ratio(i)=1;
         else
           ratio(i)=f(i)/A(i);
         end
       end
    else
       n=length(f);
       ratio=f./A;
       a1=[];
    end 
end
[a3,b1]=sort(ratio,'descend');  %将其重要性排序
total=0;
y=zeros(1,n);
for i=1:n
   if total+A(b1(i))<=b
       y(b1(i))=1;
       total=total+A(b1(i));
   else
       y(b1(i))=(b-total)/A(b1(i));
       break;
   end
end
if  ~isempty(a1)
    y(I(a1,1))=0;
end
val=f*y';     
    
    
    

