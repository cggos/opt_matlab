function y=isin_TSP(y,NC)   %将重复的城市划去，补上缺的城市
if nargin==1
    NC=length(y);
end
lost=[];
for i=1:NC
    if isin(i,y)==0           %缺的城市
      lost=[lost i];
    end
end
m=find(y>NC);    %超出范围
if ~isempty(m)
   m=[m find(y<=0)];
else
   m=find(y<=0);
end
if ~isempty(m)
    for i=1:length(m)
        y(m(i))=lost(i);   
    end
    lost=redu(lost,1:length(m),'c');
end
if ~isempty(lost)
     for i=1:NC
       a=find(y==i);
       if length(a)>1
           for j=1:length(a)-1
               y(a(j))=lost(j);
           end
           lost=redu(lost,1:length(a)-1,'c');
       end
       if isempty(lost)
           break
       end
     end 
end
if NC<length(y)
   for i=1:NC
      m=find(y==i);
      if length(m)>1
          y=redu(y,m(2:end),'c');
      end
   end
end
       