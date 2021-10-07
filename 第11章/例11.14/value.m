function fval=value(route,d)%used for reckoning the goal value of the selected traveling route
n=length(d);
fval=0;
for i=1:n
    if i==n
        fval=fval+d(route(i),route(1));
    else
        fval=fval+d(route(i),route(i+1));
    end
end