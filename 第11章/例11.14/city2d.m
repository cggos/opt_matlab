function [num,d]=city2d(city)   
[r,c]=size(city);
if r==2
    num=c;
    d=squareform(pdist(city));
elseif c==2
    num=r;
    d=squareform(pdist(city));
else
    num=r;
    d=city;
end