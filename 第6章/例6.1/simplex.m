function [xsol, fval, lambda, exitflag, output, basicVarIdx, nonbasicVarIdx, delrows] = simplex(c, Ale, ble, ...
                     Aeq, beq, lb, ub, options, defaultopt, computeLambda)
%SIMPLEX Simplex method for general linear programming problem.
%               min  c' * x 
%        subject to  Ale * x <= ble; 
%                    Aeq * x  = beq; 
%                    lb <= x <= ub.
%
%   XSOL = SIMPLEX(c,A,b,Aeq,beq) solves the problem above while 
%   additionally satisfying the equality constraints Aeq * xsol = beq.
%
%   XSOL = SIMPLEX(c,A,b,Aeq,beq,LB,UB) defines a set of lower and upper 
%   bounds on the design variables, XSOL, so that the solution XSOL is always in the
%   range LB <= XSOL <= UB.  If LB is [] or if any of the entries in LB are
%   -Inf, that is taken to mean the corresponding variable is not bounded below.
%   If UB is [] or if any of the entries in UB are Inf, that is taken to mean
%   the corresponding variable is not bounded above.
%
%   XSOL = SIMPLEX(c,A,b,Aeq,beq,LB,UB,OPTIONS) minimizes with the default 
%   optimization parameters replaced by values in the structure OPTIONS, an
%   argument created with the OPTIMSET function.  See OPTIMSET for details.
%   Used options are Display, TolFun, LargeScale, MaxIter. Use OPTIONS = [] as
%   a place holder if no options are set.
%
%   [XSOL,FVAL] = SIMPLEX(c,A,b) returns the value of the objective
%   function at XSOL: FVAL = c' * XSOL.
%
%   [XSOL,FVAL,LAMBDA] = SIMPLEX(c,A,b) returns the set of Lagrange multipliers,
%   LAMBDA, at the solution, XSOL.
%   LAMBDA.ineqlin are the Lagrange multipliers for the inequalities.
%   LAMBDA.eqlin are the Lagrange multipliers for the equalities.
%   LAMBDA.lb the Lagrange multipliers for the lower bounds.
%   LAMBDA.ub are the Lagrange multipliers for the upper bounds.
%
%   [XSOL,FVAL,LAMBDA,EXITFLAG] = SIMPLEX(c,A,b) describes the exit conditions:
%   If EXITFLAG is:
%      > 0 then SIMPLEX converged with a solution XSOL.
%      = 0 then SIMPLEX did not converge within the maximum number of iterations.
%      < 0 then the problem was infeasible or unbounded.
%
%   [XSOL,FVAL,LAMBDA,EXITFLAG,OUTPUT] = SIMPLEX(c,A,b) returns a structure:
%     OUTPUT.iterations: number of iterations.
%     OUTPUT.cgiterations: [].

%   Copyright 1990-2015 The MathWorks, Inc.

% Handle missing arguments
if nargin < 10  
    computeLambda = false;
end

if nargin < 8, options = []; 
   if nargin < 7, ub = []; 
      if nargin < 6, lb = []; 
         if nargin < 5, beq = [];
            if nargin < 4, Aeq = [];
               end, end, end, end, end

if nargout > 5
   restoreBasis = true;
else 
   restoreBasis = false;
end

% Set the default output fields 
output.iterations = 0;
output.algorithm = 'simplex';
output.cgiterations = [];
exitflag = 1;
basicVarIdx = []; % indices of basic variables
nonbasicVarIdx = []; % indices of nonbasic variables
delrows = []; % logical vector indicating deleted rows

% Options setup: tol and verbosity
tol = optimget(options,'TolFun',defaultopt,'fast');
switch optimget(options,'Display',defaultopt,'fast')
    case {'none','off'}
        verbosity = 0;  % error only
    case 'iter'
        verbosity = 2;  % iter and final
    case 'final'
        verbosity = 1;  
    case 'testing'
        verbosity = 4;  % all
    otherwise
        verbosity = 1;  % final message
end

% Check the input data
[c, Ale, ble, Aeq, beq, lb, ub, exitflagDatacheck] = datacheck(c, Ale, ble, Aeq, beq, lb, ub, verbosity);
n_orig = length(c);
nslacks = length(ble);
meq = length(beq);

maxiter = optimget(options,'MaxIter',defaultopt,'fast');
if ischar(maxiter)
  if isequal(lower(maxiter),'100*numberofvariables')
    maxiter = 100*n_orig;
  else
    error(message('optim:simplex:InvalidMaxIter'))
  end
end

% Final output messages:
presolveTermMsg = sprintf('Optimization terminated during preprocessing phase.');
presolveOptMsg = sprintf('Optimization terminated during preprocessing phase.');
presolveInfsMsg = sprintf('Exiting: The constraints are overly stringent; no feasible point exists.');
unboundedMsg = sprintf('Exiting: The problem is unbounded; the constraints are not restrictive enough.');
excdMaxiterPhase1Msg = sprintf(['Exiting: Maximum number of iterations exceeded in Phase 1;' ...
                          '\n' '         increase options.MaxIterations.']);
phaseoneInfsMsg = sprintf('Exiting: The constraints are overly stringent; no feasible starting point found.');
excdMaxiterPhase2Msg = sprintf('Exiting: Maximum number of iterations exceeded; increase options.MaxIterations.');
optMsg = sprintf('Optimization terminated.');

% Preprocessing the problem.
% First transform into the standard LP form with equality constraints and bounds by adding slack variables. 
% The variable ndel represents the number of original variables deleted from the presolve procedure.
[c0, cc, AA, bb, lbb, ubb, ndel, resolveSol, lambda, lamindx, exitflagp, MarkVars, MarkConstr] = ...
    simplexpresolve(c, Ale, ble, Aeq, beq, lb, ub, verbosity, computeLambda, restoreBasis);

if exitflagp < 0  
    exitflag =  exitflagp;
    if exitflag == -1         % This value will be remapped to -2 in linprog.m
      msg = sprintf([presolveTermMsg '\n' presolveInfsMsg]);
    elseif exitflag == -2     % This value will be remapped to -3 in linprog.m
      msg = sprintf([presolveTermMsg '\n' unboundedMsg]);
    else
      error(message('optim:simplex:WrongStatus'))
    end
    if verbosity > 0
      disp(msg)
    end
    % At presolve stage, nothing to return to output yet.
    xsol = [];
    fval = []; 
    lambda.ineqlin = [];
    lambda.eqlin   = [];
    lambda.lower   = [];
    lambda.upper   = [];
    output.message = msg;
    return;
elseif exitflagp == 0 
    % Constraint matrix is empty, all variables are fixed, singleton, etc -
    % no more variables left to solve for (presolve solves the problem
    % already); call postsolve to restore the solution and the lambda to
    % the original system.
    xs = zeros(0,1);
    fv = 0.0;
    dualvars = struct('y', [], 'z', [], 'w', []);
    caseflag1 = 1; % Empty constraint matrix
    [xsol, fval, lambda, exitflagps, basicVarIdx, nonbasicVarIdx, delrows] = simplexpostsolve(xs, fv, c0, c, resolveSol, verbosity, ...
        lambda, lamindx, dualvars, Ale, Aeq, computeLambda, caseflag1, MarkVars, basicVarIdx, nonbasicVarIdx, restoreBasis, MarkConstr, []);
    exitflag = exitflagps;

    if (exitflag == 1) 
      msg = sprintf(['Solution determined by the constraints.\n' presolveOptMsg]);
    elseif (exitflag == -1) % This value will be remapped to -2 in linprog.m
      msg = sprintf([presolveTermMsg '\n' presolveInfsMsg]);
    else
      error(message('optim:simplex:WrongStatus'));
    end
    % Display the final statement
    if verbosity > 0
      disp(msg)
    end
    output.message = msg;
    return;
end

%Step 3. phase 1 to find the first feasible basis by solving an auxiliary piecewise linear problem 
[c1, A, b, lbs, ubs, basicVarIdx, nonbasicVarIdx, x1opt, exitflagPhase1, lamindx, delrows1] = ...
    simplexphaseone(cc, AA, bb, lbb, ubb, n_orig, ndel, nslacks, maxiter, tol, verbosity, lamindx, computeLambda); 
exitflag = exitflagPhase1;
% exitflagPhase1 = 1, feasible; 
% exitflagPhase1 = 0, exceed the maximum iteration; 
% exitflagPhase1 = -1, infeasible; 
% exitflagPhase1 = -2, unbounded.

% Note that empty objective coefficients is recasted as zeros(n_orig, 1) in linprog.
if (exitflagPhase1 == 1) && all(c == 0)
    computeLambda = false;
    if (verbosity >=2)
        disp(getString(message('optim:simplex:SkipPhase2')));
    end
end

% Step 4. phase 2 is the main iterations to achieve optimality 
if exitflagPhase1 == 1 && any(c ~= 0)
    [xs, fv, dualvars, exitflagPhase2, niters, basicVarIdx, nonbasicVarIdx] = simplexphasetwo(c1, A, b, lbs, ubs, basicVarIdx, nonbasicVarIdx, x1opt, maxiter, tol, verbosity, computeLambda);
    output.iterations = niters;
    exitflag = exitflagPhase2;
    % exitflagPhase2 = 1, optimal; exitflagPhase2 =0, maxiter is exceeded in Phase 2;
    % exitflagPhase2 = -1, infeasible; exitflagPhase2 = -2, unbounded.
else
    xs = x1opt;
    fv = c1' * x1opt;
    dualvars = struct('y', [], 'z', [], 'w', []);
end 

% Step 5. postsolve restores the solution and the lambda (if applicable) to the original system                                          
caseflag2 = 2; % Nonempty constraint matrix
if exitflag <= 0
    computeLambda = false;
    restoreBasis  = false;
end
[xsol, fval, lambda, exitflagps, basicVarIdx, nonbasicVarIdx, delrows] = simplexpostsolve(xs, fv, c0, c, resolveSol, verbosity, ...
    lambda, lamindx, dualvars, Ale, Aeq, computeLambda, caseflag2, MarkVars, basicVarIdx, nonbasicVarIdx, restoreBasis, MarkConstr, delrows1);
if exitflag == 1
    exitflag = exitflagps;
end

% Make sure output vector is full.
xsol = full(xsol);

% The values -1 and -2 of exitflag will be remapped to -2 and -3, 
% respectively, in linprog.m
if exitflag == 1
  msg = optMsg;
elseif exitflagPhase1 == 0 % Note exitflag == 0.
  msg = excdMaxiterPhase1Msg;
elseif exitflag == 0
  msg = excdMaxiterPhase2Msg;
elseif exitflagPhase1 == -1 % Note exitflag == -1.
  msg = phaseoneInfsMsg;
elseif exitflag == -1
  msg = presolveInfsMsg;
elseif exitflag == -2
  msg = unboundedMsg;
else
  error(message('optim:simplex:WrongStatus'));
end
% Display the final statement
if verbosity > 0
  disp(msg)
end
output.message = msg;

%--------------------------------------------------------------------------
%-----------------------SIMPLEX Subfunctions----------------------------
%--------------------------------------------------------------------------
function [c, Ale, ble, Aeq, beq, lb, ub, exitflagDatacheck] = ...
    datacheck(c, Ale, ble, Aeq, beq, lb, ub, verbosity)
% DATACHECK Check input data on the consistency of dimensions and boundary conditions.

exitflagDatacheck = 1;
% Check input data
if (nargin < 3), ble = [];
    if (nargin < 2), Ale = [];
        if (nargin < 1) 
            error(message('optim:simplex:NotEnoughInputs'));
        end, end, end

[mle, nle] = size(Ale);
[meq, neq] = size(Aeq);
nvars = max( [length(c), length(ub), length(lb), neq, nle] );

if isempty(ub)
    ub = inf(nvars, 1);
    if (verbosity >= 4)
        disp('  Assuming all variables are unbounded above.');
    end
end

if isempty(lb)
    lb = -inf(nvars, 1);
    if (verbosity >= 4)
        disp('  Assuming all variables are unbounded below.');
    end
end  

if isempty(c), c=zeros(nvars,1); end
if isempty(Ale), Ale=zeros(0,nvars); end
if isempty(ble), ble=zeros(0,1); end       
if isempty(Aeq), Aeq=zeros(0,nvars); end
if isempty(beq), beq=zeros(0,1); end       

c = c(:);
lb = lb(:);
ub = ub(:);
ble = ble(:);
beq = beq(:);

if (size(c,2) ~= 1)  
    error(message('optim:simplex:InvalidFirstInput'))
end

if any(isnan([c; lb; ub])) || any(isnan([ble;beq])) || any(any(isnan([Ale; Aeq])))
    exitflagDatacheck = -1; 
    error(message('optim:simplex:NaNInput'));
end

if  any(any(isinf([c'; Ale; Aeq]))) || any(isinf([ble;beq]))
    exitflagDatacheck = -1;
    error(message('optim:simplex:InfInput'));
end

%--------------------------------------------------------------------------
%-------------END of SIMPLEX and its subfunctions-----------------------
%--------------------------------------------------------------------------
