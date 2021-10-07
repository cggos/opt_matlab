function [cc, A, b, lbs, ubs, basicVarIdx, nonbasicVarIdx, x1opt, exitflag, lamindx, delrows] = ...
simplexphaseone(ct, At, bt, lbt, ubt, n_orig, ndel, nslacks, maxiter, tol, verbosity, lamindx, computeLambda)
%SIMPLEXPHASEONE Find the first feasible basis solution.
%   The feasible basis solution is found by constructing and solving 
%   an auxiliary piecewise LP problem.

%   Reference:  Robert E. Bixby (1992), "Implementing the Simplex Method: 
%   The Initial Basis", ORSA Journal on Computing Vol. 4, No. 3.

%   Copyright 1990-2010 The MathWorks, Inc.

% Initialization
tol2 = 1.0e-2 * tol;
A = At;
b = bt;
lbs = lbt;
ubs = ubt;
lb = lbt;
ub = ubt; 
cc = ct;
c = cc;
c = c(:);
[m, nn] = size(A);  
delrows = [];

%----Order the original variables x(1:n) according to the preference-----

VarSetlf = isfinite(lb) & (ub == Inf);
VarSetuf = isfinite(ub) & (lb == -Inf);
VarSetf  = isfinite(lb) & isfinite(ub);

% Define the penalty vector p and sort the columns accordingly
p = zeros(nn,1);
p(VarSetlf) = lb(VarSetlf);
p(VarSetuf) = - ub(VarSetuf);
p(VarSetf)  = lb(VarSetf) - ub(VarSetf);
gamma = max(abs(c));
if gamma ~= 0 
    cmax = gamma * 1000;
else
    cmax = 1;
end
p = p + c/cmax;
[p, idxorder] = sort(p);

% Step 1 Initialization of the first basis
% For the rows corresponding to inequalities, set I(i) and r(i) to 1.
% The nbasicVarIdx is the current number of basis variables(columns) we have chosen.
% Leftover slack variables are chosen to be basis variables.
basicVarIdx =  ( (n_orig-ndel+1): nn )';
nbasicVarIdx = nnz(basicVarIdx);
I = ones(m,1); 
I( (nbasicVarIdx +1):m) = 0;                     
r = I;             
pvt = inf(m,1);

nonbasicVarIdx = true(nn,1);
nonbasicVarIdx(basicVarIdx) = false ; % Assign basis vectors to false in nonbasis

% Step 2 Iteration on each column
for j = idxorder'     
    ii = find(r == 0);
    [alpha, maxi] = max( abs(A(ii, j)) );
    maxidx = ii(maxi); 
    if ~isempty(alpha) && (alpha >= 0.99) && ( nbasicVarIdx + 1 <= m)
        nbasicVarIdx = nbasicVarIdx + 1;
        basicVarIdx( nbasicVarIdx, 1 ) = j; % nonbasicVarIdx changes.
        nonbasicVarIdx(j,1) = false;     % Mark column x_j as basis variable in nonbasis.
        
        I(maxidx) = 1;
        pvt(maxidx) = alpha;
        ii = ( abs(A(:, j)) ~= 0 );
        r(ii) = r(ii) + 1;
        break;
    end
    
    if any( abs(A(:, j )) > 0.01*pvt )
        break;
    else
        ii= find(I == 0);
        [alpha, maxi] = max( abs(A(ii, j)) );
        maxidx = ii(maxi);
        if ( alpha == 0 )
            break;
        elseif (nbasicVarIdx +1 <= m)
            nbasicVarIdx = nbasicVarIdx + 1;
            basicVarIdx( nbasicVarIdx, 1 ) = j;
            nonbasicVarIdx(j, 1) = false;
            
            I(maxidx) = 1;
            pvt(maxidx) = alpha;
            ii = ( abs(A(:, j)) ~= 0 );
            r(ii) = r(ii) + 1;
        end
    end
end % for j = idxorder'

% Step 3. Add artificial variables when necessary.
% the number of rows uncovered = the number of artificial variables to add
ii = sum(I); 
nz = m - ii;  %  Here m is the current number of rows.
if nz >= 1
    ii = find(I == 0);
    A = [A sparse(ii, 1:nz, ones(nz,1), m, nz)];
    basicVarIdx((nbasicVarIdx+1):(nbasicVarIdx +nz), 1)= ((nn+1) :(nn+nz))';
end

B = A(:, basicVarIdx);
Nfix = nonbasicVarIdx & (lb == ub); 
Nb   = nonbasicVarIdx & (lb ~= ub) & (isfinite(lb) | isfinite(ub) );
Nlb  = Nb & ( abs(lb) <= abs(ub) );   % Nlb is a logical vector nonbasicVarIdx(Nlb) 
Nub  = Nb & ( abs(lb) >  abs(ub) );   % Nub is a logical vector nonbasicVarIdx(Nub)

% Obtain the basic solution (may not be feasible) to the auxiliary problem
y = zeros(nn, 1);
y(Nfix) = lb(Nfix);
y(Nlb)  = lb(Nlb);
y(Nub)  = ub(Nub);
y = [y; zeros(nz,1)];
An  = A(:, 1:nn);
if any(nonbasicVarIdx)
    rhs = b - An(:, Nlb)*y(Nlb) - An(:, Nub)*y(Nub) - An(:, Nfix)*y(Nfix);
else
    rhs = b;
end

% Disable the warnings about conditioning for singular and
% nearly singular matrices
warningstate1 = warning('off', 'MATLAB:nearlySingularMatrix');
warningstate2 = warning('off', 'MATLAB:singularMatrix');

y(basicVarIdx)= B\rhs;      

% Restore the warning states to their original settings
warning(warningstate1)
warning(warningstate2)

ylb = [lb;  zeros(nz,1)];
yub = [ub;  zeros(nz,1)];

% Set up the objective of the auxiliary problem based on the basis solution assumptions: 
% p is the penalty coefficient for the auxiliary problem with constraint matrix Ay = b 
% y is the solution to the initial basis for the auxiliary problem

p = zeros(nn +nz, 1);
pconst = 0; 
llb = y < ylb;
uub = y > yub;
if any(llb)
    p(llb) = -1;
    pconst = pconst + sum(ylb(llb));
end
if any(uub)
    p(uub) = 1; 
    pconst = pconst - sum(yub(uub)); 
end

exitflag = 1; % Some of the values of exitflag set in this file will be remapped
              % in linprog.m
ii = (1:nn)';                     
nonbasicVarIdx = setdiff(ii(nonbasicVarIdx), basicVarIdx); % Return the value of nonbasicVarIdx instead of logical type

% Call the iteration procedure SIMPLEXPIECEWISE to solve the auxiliary problem if necessary
if all(p == 0)
    x1opt = [y(basicVarIdx); y(nonbasicVarIdx)]; 
    exitflag = 1;
    if verbosity >= 2
        disp(getString(message('optim:simplexphaseone:FeasStartPtSkipPhase1')));
    end
else
    if (verbosity >= 4)
        disp( sprintf('  Call the simplexpiecewise procedure to solve the auxiliary problem.') );    
        % For checking purpose
        checkin1 = any( abs(A*y - b) > tol2 );
        if (checkin1)
            error(message('optim:simplexphaseone:WrongSetup1'));
        end
    end
    
    [x1opt, pval, basicVarIdx, nonbasicVarIdx, exitflag, niters] = simplexpiecewise(n_orig, -p, A, b, ylb, yub, basicVarIdx, nonbasicVarIdx, y, maxiter, tol, verbosity);
    
    % Assumption: simplexpiecewise returns x1opt = [x_B, x_N].
    boundcheck = all( x1opt <= yub([basicVarIdx; nonbasicVarIdx]) + tol2 ) && all(x1opt >= ylb([basicVarIdx; nonbasicVarIdx]) - tol2);
    constrcheck = all( abs(A(:, [basicVarIdx; nonbasicVarIdx])* x1opt - b) <= tol2); 
    
    if ( abs(pval) <= 1.e-10 ) && boundcheck && constrcheck
        exitflag = 1;
        if verbosity >= 4
            disp( sprintf('  Find the first feasible basis.') );
        end
    elseif (exitflag ~= 0)
        exitflag = -1; % This value will be remapped to -2 in linprog.m
        if (verbosity >= 4)
            disp( sprintf('The problem is infeasible, detected by the auxiliary problem.') );
        end
    end
end

% Order x1opt according to the column of A before function return
if ((exitflag == 1) && (nz == 0)) || (exitflag <= 0)
    [tmp, order] = sort( [basicVarIdx; nonbasicVarIdx] );
    x1opt = x1opt(order);
    x1opt = x1opt(1: length(cc));
    return;     
end

%------The procedure of getting rid of artificial basic variables-------- 
% Both a feasible solution and artificial variables exist.
if (exitflag == 1) &&  (nz > 0) 
    
    delcols  = false(1, size(A,2) );
    mB       = size( basicVarIdx, 1 );
    delB     = false(1, mB);
    artBVars = find( basicVarIdx > nn );
    
    if ( ~isempty(artBVars) )
        delrows  = false(mB, 1);
        % Disable the warnings about conditioning for singular and
        % nearly singular matrices
        warningstate1 = warning('off', 'MATLAB:nearlySingularMatrix');
        warningstate2 = warning('off', 'MATLAB:singularMatrix');
        for k = 1:length(artBVars)
            % Pick the one to leave the Basis.
            % e is the leaving row of the identity matrix
            % Solve the linear system r'*B = e(lv), where
            % r is a row vector to be determined.
            % compute rN = r*N 
            lv    = artBVars(k);
            e     = zeros(1, size(A,1));
            e(lv) = 1;
            r     = e / A(:, basicVarIdx); 
            ind   = find(nonbasicVarIdx <= nn);    
            rN    = r * A(:, nonbasicVarIdx(ind));
            
            candEnt = ind( abs(rN) > tol2 ); 
            if  isempty(candEnt)             
                % Then the corresponding row of constraint matrix A is redundant.
                if ( abs(r*b) > tol2) && (verbosity >= 4) % For double-check purpose 
                    % Restore the warning states to their original settings
                    warning(warningstate1)
                    warning(warningstate2)
                    error(message('optim:simplexphaseone:Inconsistency'));
                end
                delrows( A(:, basicVarIdx(lv)) ~= 0 ) = true;
                delB(lv) = true;
                delcols(basicVarIdx(lv)) = true;
            else                             
                ent = candEnt(1); % Pick the first candidate that is not artificial
                
                % Replace the leaving artificial variable by the entering variable
                swaptmp      = basicVarIdx(lv);
                basicVarIdx(lv)    = nonbasicVarIdx(ent);
                nonbasicVarIdx(ent)   = swaptmp;
                
                swaptmp      = x1opt(lv);
                x1opt(lv)    = x1opt(mB + ent);
                x1opt(mB+ent)= swaptmp;
            end % if ~any(candEnt)
            
        end % for k=1: length(artBvars)
        % Restore the warning states to their original settings
        warning(warningstate1)
        warning(warningstate2)
        
        % The remaining artificial basic variables after the above procedure (for loop) 
        % corresponds to redundant rows in equality constraints Ax = b.
        % Delete redundant rows in A and b and delete the remaining artificial columns.
        if any(delrows)
            A(delrows, :) = [];
            b(delrows)    = [];
            basicVarIdx(delB) = [];
            if computeLambda
                Aindex = lamindx.Aindex;
                leftrows = ~delrows;
                lamindx.Aindex = Aindex(leftrows);
                lamindx.yindex = leftrows; % to signal the deleted rows from rowscl
            end
            if verbosity >= 4
                disp( sprintf('  Delete %d redundant rows in the constraint matrix in SIMPLEXPHASEONE.',...
                    nnz(delrows) ) );
            end
        end
        
    end % if ~isempty(artBVars)
    
    % Delete the artificial nonbasic variables
    artNVars = ( nonbasicVarIdx > nn );
    delcols( nonbasicVarIdx(artNVars) ) = true;
    A(:, delcols ) = [];
    
    nonbasicVarIdx( artNVars ) = [];
    x1opt( [delB'; artNVars] ) = [];
    
    % Order the output of x1opt corresponding to the original problem
    [tmp, order] = sort( [basicVarIdx; nonbasicVarIdx] );
    x1opt = x1opt(order);
    
    if verbosity >= 4
        % check the feasibility of  Ax = b, lb <= x1opt <= ub 
        checkin2 = any( abs(A*x1opt - b) > tol2 );
        if (checkin2)
            error(message('optim:simplexphaseone:WrongSetup2'));
        end
        
        ubcheck = max ( x1opt - ubs );
        lbcheck = max ( lbs - x1opt );
        if (ubcheck > tol2) || (lbcheck > tol2)
            error(message('optim:simplexphaseone:InvalidBndryCond', sprintf( '%8.2e', max( ubcheck, lbcheck ) )));
        end
    end
    
    if verbosity >= 4
        % Check the output data 
        [mA, nA] = size(A);
        if mA ~= length(b)
            error(message('optim:simplexphaseone:SizeABMismatch'));
        end
        
        if nA ~= length(cc)
            error(message('optim:simplexphaseone:SizeACcMismatch'));
        end
        
        if nA ~= length(lbs)
            error(message('optim:simplexphaseone:SizeALbsMismatch'));
        end
        
        if nA ~= length(ubs)
            error(message('optim:simplexphaseone:SizeAUbsMismatch'));
        end
        
        if nA ~= length(x1opt)
            error(message('optim:simplexphaseone:SizeAX1optMismatch'));
        end
        
        if nA ~= ( length(basicVarIdx)+length(nonbasicVarIdx) )
            error(message('optim:simplexphaseone:SizeAVarIdxMismatch'));
        end 
    end % if verbosity >= 4
    
end % if (exitflag == 1) &&  (nz > 0) end of the procedure of deleting artificial variables
