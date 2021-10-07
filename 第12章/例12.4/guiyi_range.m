function y=guiyi_range(x,x1)  %归一化到某一区间
[r,c]=size(x);
if r==1
    x=x';
end
if ~iscell(x)
  y_min=min(x);
  y_max=max(x);
  R=y_max-y_min;
  for i=1:r
    for j=1:c
        y(i,j)=x1(1)+(x1(2)-x1(1))*(x(i,j)-y_min(j))/R(j);
    end
  end
else
   y_min=min(cell2mat(x));
   y_max=max(cell2mat(x));
   R=y_max-y_min;
   for i=1:r
     for j=1:c
        y(i,j)={x1(1)+(x1(2)-x1(1))*(cell2mat(x(i,j))-y_min(j))/R(j)};
     end
   end
end
