function H=myhessian1(fun,x_syms)    %求海赛矩阵,x_syms为变量
if isempty(fun)
    H=[];
else
  n=size(fun,1);
  num=length(x_syms);
  H=[];
  for k=1:n
     f=fun(k,:);
     for i=1:num
         d=diff(f,x_syms(i));
         for j=1:num
            y(i,j)=diff(d,x_syms(j));
         end
     end
     H=[H;y];
  end
end
