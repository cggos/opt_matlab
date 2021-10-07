function r = simplify(s,varargin)
%SIMPLIFY Symbolic simplification.
%   SIMPLIFY(S) simplifies each element of the symbolic matrix S.
%   
%   SIMPLIFY(S,N) or, equivalently, SIMPLIFY(S,'Steps',N),
%   searches for a simplification in N steps. The default value of N is 1.
%   
%   SIMPLIFY(S,'Seconds',T) aborts the search for a simpler version
%   at the next end of a simplification step after T seconds. The results
%   will depend on the speed and system load of your computer and may
%   not be reproducible.
%   
%   SIMPLIFY(S,'IgnoreAnalyticConstraints',VAL) controls the level of 
%   mathematical rigor to use on the analytical constraints while simplifying 
%   (branch cuts, division by zero, etc). The options for VAL are TRUE or 
%   FALSE. Specify TRUE to relax the level of mathematical rigor
%   in the rewriting process. The default is FALSE.
%   
%   SIMPLIFY(S, 'Criterion', 'preferReal') 
%   discourages simplify from returning results 
%   containing complex numbers.
%
%   Examples: 
%      simplify(sin(x)^2 + cos(x)^2)  returns  1.
%      simplify(exp(c*log(sqrt(alpha+beta))))  returns  (alpha + beta)^(c/2).
%      simplify(sqrt(x^2))  returns  sqrt(x^2),
%      simplify(sqrt(x^2),'IgnoreAnalyticConstraints',true)  returns  x.
%
%   See also SYM/FACTOR, SYM/EXPAND, SYM/COLLECT.

%   Copyright 1993-2014 The MathWorks, Inc.

if builtin('numel',s) ~= 1,  s = normalizesym(s);  end

posintvalidator = @(value) isnumeric(value) && 0 <= value && round(value) == value;

p = inputParser;
p.addOptional('Steps', 1, posintvalidator);
p.addParameter('Seconds', sym(inf), @(value) isnumeric(value) && 0 <= value);
p.addParameter('Criterion', 'default');
p.addParameter('IgnoreAnalyticConstraints', false);
p.parse(varargin{:});

steps = sym(p.Results.Steps);
seconds = sym(p.Results.Seconds);
criterion = validatestring(p.Results.Criterion, {'default', 'preferReal'});
if strcmp(criterion, 'preferReal')
    criterion = evalin(symengine, 'Simplify::preferReal');
else
    criterion = evalin(symengine, 'Simplify::defaultValuation');
end    
if sym.checkIgnoreAnalyticConstraintsValue(p.Results.IgnoreAnalyticConstraints) 
   IAC = 'TRUE';
else
   IAC = 'FALSE';
end   

rSym = feval(symengine, 'simplify', s,...
            evalin(symengine, ['IgnoreAnalyticConstraints =', IAC]),...
            sym('Steps') == steps,...
            sym('Seconds') == seconds,...
            sym('Criterion') == criterion,...
            evalin(symengine, 'OutputType = "Best"') ...
            );
r = privResolveOutput(rSym, s);
