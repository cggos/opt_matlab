function [y,dy]=opfun2(varargin)
if isnumeric(varargin{1})
   x=varargin{1};
   y=3*x^2-2*tan(x);
   dy=[];
else
   y=eval(subs(varargin{1},varargin{2},varargin{3}));
   y1=mydiff(varargin{1},varargin{2});
   dy=eval(subs(y1,varargin{2},varargin{3}));
end