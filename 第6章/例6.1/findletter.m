function [y,num]=findletter(str,type,symbol)    %从字符中找到字母
if nargin==2
   symbol=',';
end
if isempty(find(char(str)==symbol))
    y=str;
    num=size(y,2);
else
   if type==1
      a=find(char(str)==symbol);
      y=redu(str,a,'c');
      num=size(y,2);
   elseif type==2   %有空格时
      if ~strcmp(symbol,' ')
         b=findstr(char(str),' ');
         if ~isempty(b)
           str=redu(str,b,'c');
         end
      end
      a=find(char(str)==symbol);
      n=length(a);
      for i=1:n
         if i==1
            y{i}=char(str(1:a(i)-1));
         else
            y{i}=char(str(a(i-1)+1:a(i)-1));
         end
      end
      if a(end)<length(str)
         y{n+1}=char(str(a(end)+1:end));
         num=n+1;
      else
         num=n;
      end
   end 
end
    
