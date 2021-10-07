function y=TSPplot(city,route)
n=length(route);
plot(city(:,1),city(:,2),'o');
hold on
for i=1:n+1
    if i<=n
       x(i)=city(route(i),1);
       y(i)=city(route(i),2);
    else
        x(i)=city(route(1),1);
        y(i)=city(route(1),2);
    end
end
y=plot(x,y,'o-');
if n<=30
   gname;
end