function y=isin_TSP(y)   %���ظ��ĳ��л�ȥ������ȱ�ĳ���
NC=length(y);
lost=[];
for i=1:NC
    if isin(i,y)==0           %ȱ�ĳ���
      lost=[lost i];
    end
end
m=find(y>NC);    %������Χ
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

       