function y=mydec(x)     %��x����С������
n=length(x);
for i=1:n
  y(i)=x(i)-floor(x(i));
end