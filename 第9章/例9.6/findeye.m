function y=findeye(x)     %找矩阵x中的单位阵
[r,c1]=size(x);
m=[];
for i=1:c1
    a=find(x(:,i)==1);
    if ~isempty(a)&&length(a)==1
        m=[m i];
    end
end
if isempty(m)
    y=[];
else
   y2=[];
   for i=1:length(m)
      a=find(x(:,m(i))==1);
      b=redu(x(:,m(i)),a,'r');
      y1=find(abs(b)<0.0001);
      if length(y1)==r-1
        y2=[y2;m(i) a];
      end  
   end
   if ~isempty(y2)
      [a,b]=sort(y2(:,2));
      y2=y2(b,:); 
      y=y2;
   else
     y=[];
   end
end


    

        
    