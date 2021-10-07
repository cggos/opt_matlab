function y=myjacobian1(fun,x_syms)   %�󸴺����ĵ���,x_symsӦΪ����
if nargin==1
    syms x
    x_syms=sym(x);
end
if ~isempty(fun)
  r=size(fun,1);
  n=length(x_syms);
  for i=1:r
     for j=1:n
         y(i,j)=jacobian(fun(i,:),x_syms(j));
     end
  end
else
    y=[];
end