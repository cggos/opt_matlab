function y=sum1(x,ndim)
if nargin<2
    ndim=1;
end
[r,c]=size(x);
if r==1
   y=0;
   for i=1:c
      y=y+x(i)^2;
   end
elseif ndim==1  %按列计算
     y=zeros(c,1);
     for j=1:c
        for i=1:r
           y(j)=y(j)+x(i,j)^2;
        end
     end
elseif ndim==2    %按行计算
    y=zeros(r,1);
    for i=1:r
       for j=1:c
          y(i)=y(i)+x(i,j)^2;
       end
    end
end