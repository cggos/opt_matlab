 function y=my_delrow(y,x)
 for i=1:length(x)   
    y(find(y==x(i)))=[];
 end