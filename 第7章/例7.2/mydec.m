function y=mydec(x)     %求x的正小数部分
n=length(x);
for i=1:n
  y(i)=x(i)-floor(x(i));
end