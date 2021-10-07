function y=findzeros(x,s)   %找矩阵中为's'行或列
if nargin==1
    s=0;
end
if ~iscell(x)
  if s==0
     a=find(abs(x)<=1e-4);
  else
     a=find(x==s);
  end
  if ~isempty(a)&&length(a)==length(x)
     y='all';
  elseif ~isempty(a)
     y=a;
  else
     y=[];
  end
else
  n=length(x);y=[];
  for i=1:n
     if s==0
       a=find(abs(x{i})<=1e-4); 
     else
       a=find(x{i}==s);
     end
     if ~isempty(a)
       y=[y i];
     end
  end
end
      
 