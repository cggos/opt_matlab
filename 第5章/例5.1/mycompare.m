function [y1,y2]=mycompare(x1,x2)  %x1Îª³¤¾ØÕó,ÇóÒì¾ØÕó
 if length(x1)>=length(x2)
    a=ismember(x1,x2);
    b=x1;
 else
    a=ismember(x2,x1);
    b=x2;
 end
 num=length(a);
 y1=[];y2=[];
 for m=1:num
    if a(m)==0
        y1=[y1 b(m)];
        y2=[y2 m];
    end
 end
 