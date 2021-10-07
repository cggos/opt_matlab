function C=lattice(C,r1,c1)   %r1,c1指定一个零开始找独立格子集
[r,c]=size(C);
if nargin==1
   r1=[];c1=[];   
end
if nargin==3
    C(r1,c1)=-inf;
    y1=find(C(r1,:)==0);
    C(r1,y1)=inf;
    y1=find(C(:,c1)==0);
    C(y1,c1)=inf;
end
[y1,y2]=find(C==0);
y3=findcount(y1); %找每行最少的零
if y3(1,2)>=2 
   y4=find(C(y3(1,1),:)==0);
   C(y3(1),y4(1))=-inf;
   C(y3(1),y4(2:end))=inf;
   y5=find(C(:,y4(1))==0);   %
   C(y5,y4(1))=inf;
end
while 1
   y1=find(C==0);
   if isempty(y1)
       break
   else
      for i=1:r     %找每行或每列单独零
          y1=find(C(i,:)==0);
          if length(y1)==1
             C(i,y1)=-inf;
             y2=find(C(:,y1)==0);
             C(y2,y1)=inf;
          end
      end
      for i=1:c
         y1=find(C(:,i)==0);
         if length(y1)==1
            C(y1,i)=-inf;
            y2=find(C(y1,:)==0);
            C(y1,y2)=inf;
         end
      end
   end
end