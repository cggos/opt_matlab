function y=mymin(fun,varargin)   %判断是否为极值
d1=diff(fun);
d2=diff(fun,2);
if length(varargin)==2
    x0=varargin{1};type=varargin{2};
    a=ff(d2,x0,type);
    for i=1:length(x0)
       if x0(i)==a
           y(i)=1;
       else
           y(i)=0;
       end
    end   
elseif length(varargin)==1
   type=varargin{1}; 
   x0=eval(solve(d1));
   y=ff(d2,x0,type);
end

function y1=ff(d2,x1,type)
for i=1:length(x1)
    a=subs(d2,x1(i));
    switch type
        case 'min'
          if a>0
            y1(i)=x1(i);
          end
        case 'max'
          if a<0
            y1(i)=x1(i);
          end
    end
end
if exist('y1','var')==0
    y1=[];
end
