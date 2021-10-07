function y=opfun9(x)
if ~isa(x,'cell')
    n=size(x,1);
    for i=1:n
       a=x(i,1);b=x(i,2);
       y(i)=4*a^2+b^2-40*a-12*b+136;
    end
else
    n=length(x);
    for i=1:n
       a=x{i}(1);b=x{i}(2);
       y(i)=4*a^2+b^2-40*a-12*b+136;
    end
end