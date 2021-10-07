function [c0curr, ccurr, Acurr, bcurr, lbcurr, ubcurr, ndel, restoreSol, ...
         lambda, lamindx, PresolExit, MarkVars, MarkConstr] = simplexpresolve(c, Aineq, bineq, ...
         Aeq, beq, lb, ub, verbosity, computeLambda, restoreBasis)
%SIMPLEXPRESOLVE Preprocessing module to reduce redundancy in constraints.
%   Input: the linear programming problem
%               min  c' * x 
%        subject to  Aineq * x <= bineq; 
%                    Aeq * x  =  beq; 
%                    lb <= x <= ub
%   
%   Output:
%        PresolExit 
%               = 1,  no infeasibility or unboundedness be detected.
%               = 0,  presolve solves the problem already.
%               = -1, primal infeasible detected from preprocessing.
%               = -2, primal unbounded or dual infeasible detected from preprocessing.

%   Copyright 1990-2014 The MathWorks, Inc.

n_orig = size(c, 1);
n = n_orig;
ndel = 0;

numRowsAineq = size(Aineq,1);
numRowsAeq = size(Aeq,1);
    
% Collect some information for the calculation of Lagrange multipliers.
if computeLambda
    % Aindex is a vector of row indices (of [Aineq;Aeq]) of remaining rows, which are stored
    % in Acurr. E.g., Aindex = [3 5 9] indicates that rows 3, 5, and 9 of [Aineq;Aeq] remain
    % in the problem, and are stored in Acurr.
    Aindex = 1:numRowsAineq+numRowsAeq; 
    % Aineqstart, Aineqend, Aeqstart, Aeqend refer to the rows in matrix [Aineq; Aeq] and
    % are static.
    Aineqstart = 1; Aineqend = 0;  % end value is zero since not in A yet
    Aeqstart = 1; Aeqend = 0;      % end value is zero since not in A yet
    lbend = length(lb);  % we don't remove the inf ones, we just ignore them
    ubend = length(ub);  % we don't remove the inf ones, we just ignore them
    lbindex = 1:lbend;
    ubindex = 1:ubend;
    lambda.ineqlin = zeros(numRowsAineq,1);
    lambda.eqlin = zeros(numRowsAeq,1);
    lambda.upper = zeros(n_orig,1); 
    lambda.lower = zeros(n_orig,1); 
else
    lambda.ineqlin = [];
    lambda.eqlin = [];
    lambda.upper = [];
    lambda.lower = [];
    lamindx.Aindex = [];
end
% Logical vector AindexEq has exactly one element per remaining row (i.e., as many elements as 
% rows in Acurr). AindexEq(i) = true if ith row in Acurr is an  equality, = false otherwise.
% We need this vector, both if lambda is requested by user or not, to locate column singletons
% in equality constraints.
AindexEq = [false(numRowsAineq,1); true(numRowsAeq,1)]; 


% Initialize variables for presolve procedure
if isempty(Aineq) % Equality constraints only;
    ccurr = c;
    Acurr = Aeq;
    bcurr = beq;
    lbcurr = lb;
    ubcurr = ub;
    mcurr = size(Acurr,1);
    ncurr = n_orig;
    nslacks = 0;
    if computeLambda
        Aeqstart = 1; Aeqend = mcurr;
    end
    if (verbosity >= 4)
        disp( sprintf(['  Linear Programming problem has %d equality' ...
                ' constraints on %d variables.'], mcurr, ncurr) );
    end
else
    mineq = size(Aineq,1);
    nslacks = mineq;
    meq = size(Aeq,1);
    mcurr = mineq + meq;
    ccurr = [c; sparse(mineq,1)];
    Acurr = [Aineq speye(mineq); Aeq sparse(meq,mineq)];
    bcurr = [bineq; beq];
    lbcurr = [lb; sparse(mineq,1)];
    ubcurr = [ub; inf(mineq,1)];
    if computeLambda
        Aineqstart = 1; 
        Aineqend = numRowsAineq;
        Aeqstart = Aineqend + 1; 
        Aeqend = Aeqstart + numRowsAeq - 1;
        lbindex = 1:length(lbcurr);
        ubindex = 1:length(ubcurr);
    end
    if (verbosity >= 4)
        disp( sprintf(['  Adding %d slack variables, one for each inequality' ...
                ' constraint, to the existing\n  set of %d variables,'], nslacks, n_orig) );
    end
    ncurr = n_orig + mineq;
    if (verbosity >= 4)
        disp( sprintf('  Resulting in %d equality constraints on %d variables.', ...
            mineq,ncurr) );
        if  ~isempty(Aeq)
            disp( sprintf(['  Combining the %d (formerly inequality) constraints' ...
                    ' with %d equality constraints,\n  resulting in %d equality' ...
                    ' constraints on %d variables.\n'], mineq, meq, mcurr, ncurr) );
        end
    end
end

% Initialize the constant part of the objective function c0curr to be zero
meq0 = mcurr;
c0curr      = 0.0;
PresolExit  = 1;

% Initialize the structure for information needed to restore the solution to the original system
restoreSol(1:10) = struct('case', '', ...
    'exist', false, ...
    'unslvd_idx', [], ...
    'slvd_idx', [], ...
    'slvd_val', [], ...
    'Adel', [], ...
    'bdel', []);

restoreSol(10).case = 'ZeroRows';
restoreSol(9).case = 'RowDependent';
restoreSol(8).case = 'FixedVarsbyObj';
restoreSol(7).case = 'ZeroColsGenerated';
restoreSol(6).case = 'FixedVarsbyBounds';
restoreSol(5).case = 'RowSgtons';
restoreSol(4).case = 'ForConstr';
restoreSol(3).case = 'ColSgtons';
restoreSol(2).case = 'ZeroCols';
restoreSol(1).case = 'ColScal';

% Tags to all constraints and variables
% MarkVars(i) = true, the ith variable remains in the system
% MarkVars(i) = false, the ith variable is deleted from the system
MarkConstr  = true(mcurr,1); 
MarkVars    = true(ncurr,1);

% Make the constraint matrix and objective coefficient sparse
if (~issparse(Acurr))
    Acurr = sparse(Acurr);
end
ccurr = sparse(ccurr);
nvars = n_orig + nslacks;

if computeLambda
    % We might delete variables and then also the corresponding components of
    % ub/lb (and also w/z) so we need a separate index into w/z 
    windex = true(length(ubcurr), 1);
    zindex = true(length(lbcurr), 1);
    yindex = true(size(Acurr,1),1);
end

% Perform preliminary solving, rearranging and checking for infeasibility

% Case: remove all fixed variables by equal bounds
fixed      = (lbcurr == ubcurr);
FixedExist = any(fixed);
if (FixedExist)
    unfix    = ~fixed;  
    xfix     = lbcurr(fixed);    
    bcurr    = bcurr - Acurr(:,fixed)*xfix;
    if  ( sum(fixed) == ncurr) && any(bcurr ~= 0)
        if verbosity >= 4
            disp( sprintf(['  Exiting due to infeasibility: Violate ' ...
                    'the equality constraints while all variables fixed.']) );
        end
        lambda.ineqlin = [];
        lambda.eqlin = [];
        lambda.upper = [];
        lambda.lower = [];
        lamindx.Aindex = [];
        PresolExit = -1; % Will be remapped to -2 in linprog.m
        return
    end
    
    if computeLambda
        % Record the deletion of variables with fixed bounds
        fixedlower = (ccurr >=0) & fixed;
        fixedupper = (ccurr < 0) & fixed;
        restoreSol(6).lorow = fixedlower;
        restoreSol(6).upcol = fixedupper;
        
        lbindex = lbindex(unfix);
        zindex = zindex(unfix); 
        ubindex = ubindex(unfix);
        windex = windex(unfix);
    end
    
    % Update 
    xpresol(fixed) = xfix;
    c0curr = c0curr + ccurr(fixed)' * xfix;
    MarkVars(fixed) = false; 
    
    ccurr  = ccurr(unfix);
    Acurr  = Acurr(:,unfix);
    lbcurr = lbcurr(unfix);
    ubcurr = ubcurr(unfix);
    ncurr  = full(sum(unfix));
    
    if (verbosity >= 4)
        disp( sprintf('  Assigning %i equal values of lower and upper bounds to solution.',...
            sum(full(fixed)) ) );
        if (ncurr == 0);
            disp('  All variables are fixed.');
        end
    end
    
    % Record the recovery info
    restoreSol(6).exist = true;
    restoreSol(6).unslvd_idx = unfix;
    restoreSol(6).slvd_idx = fixed;
    restoreSol(6).slvd_val = xfix;
end 

% Case: remove zero rows (for now in equalities only; here
% all inequalities have a nonzero coefficient corresponding 
% to the slack)
ZrowsExist = 0;
rnnzct = sum(spones(Acurr'), 1);     
if ~isempty(Acurr) && (any(rnnzct == 0))
    zrows = (rnnzct == 0);
    ZrowsExist = any(zrows);
    if (any(bcurr(zrows) ~= 0))
        if (verbosity >= 4)
            disp( sprintf(['  Exiting due to infeasibility: An all-zero row in the constraint' ...
                    ' matrix does not have a zero in corresponding right-hand-side entry.']) );
        end
        lambda.ineqlin = [];
        lambda.eqlin = [];
        lambda.upper = [];
        lambda.lower = [];
        lamindx.Aindex = [];
        PresolExit = -1; % Will be remapped to -2 in linprog.m
        return
    else   
        nzrows = ~zrows;
        Acurr  = Acurr(nzrows,:);
        bcurr  = bcurr(nzrows); 
        rnnzct = rnnzct(nzrows);
        if (verbosity >= 4)
            disp( sprintf('  Deleting %i all-zero rows from the constraint matrix.', nnz(zrows)) );
        end
        
        if computeLambda
            % Lagrange multiplier for zero rows are zero
            Aindex = Aindex(nzrows);
            yindex = yindex(nzrows);
        end
        AindexEq = AindexEq(nzrows);
       
        if restoreBasis
            rowidx = find(MarkConstr); % index of remaining rows
            MarkConstr(rowidx(zrows)) = false; % flag zero rows
            restoreSol(10).slvd_idx = rowidx(zrows); % restore zero row index
        end
        mcurr = size(bcurr,1);
    end
    
    restoreSol(10).exist = ZrowsExist;
end

% Case: make A structurally "full rank"
%       note the structure rank of the sparse matrix sprank(A) >= rank(A)
%       so dependent rows may still exist, may be picked up in simplexphaseone.m
sprk = sprank(Acurr');
Row_deleted = 0;
if (sprk < mcurr) && ~isempty(Acurr)
    Row_deleted = 1;
    [dmp, tmp] = dmperm(Acurr);
    irow = dmp(1:sprk);          % Note permutation of rows might occur. 
    idelrow = dmp(sprk+1:mcurr);
    Adel = Acurr(idelrow,:);
    bdel = bcurr(idelrow);
    iroword = sort(irow);
    Acurr = Acurr(iroword,:);
    bcurr = bcurr(iroword);
    if computeLambda
        Aindex = Aindex(iroword);
        yindex = yindex(iroword);
    end
    AindexEq = AindexEq(iroword);
    
    rnnzct = rnnzct(iroword);
    if (verbosity >= 4)
        disp( sprintf('  Deleting %i dependent rows from constraint matrix.', ... 
            mcurr-sprk ) );
    end
    mcurr = size(Acurr,1);
    if restoreBasis
        rowidx = find(MarkConstr);
        MarkConstr(rowidx(idelrow)) = false; % mark off dependent rows
    end
    
    % Record the information
    restoreSol(9).exist = true;  
    restoreSol(9).Adel = Adel;      
    restoreSol(9).bdel = bdel;
    restoreSol(9).slvd_idx = idelrow; % the index of deleted dependent rows
end

% Case: delete zero columns
ZrcolsExist = 0;
if isempty(Acurr)
    zrcol = 0; 
else
    zrcol = (max(abs(Acurr)) == 0)';  % max behaves differently if A only has one row:
                                      % so zero cols won't be deleted in this case.
end

if ( any(zrcol == 1) )
    ZrcolsExist = 1;
    izcol = find(zrcol);
    if ( any( (ccurr(zrcol) < 0) & isinf(ubcurr(zrcol)) ) ...
            ||  any( (ccurr(zrcol) > 0) & isinf(lbcurr(zrcol)) ) )
        if (verbosity >= 4)
            disp( sprintf(['  Exiting due to unboundedness: ', ...
                    'The problem is unbounded detected by zero columns.']) );
        end
        lambda.ineqlin = [];
        lambda.eqlin = [];
        lambda.upper = [];
        lambda.lower = [];
        lamindx.Aindex = [];
        PresolExit = -2; % Will be remapped to -3 in linprog.m
        return
    end
    xzrcol = zeros(nnz(zrcol), 1);
    fubindx = ccurr(zrcol) < 0;
    flbindx = ccurr(zrcol) > 0;
    tub = ubcurr(zrcol); 
    tlb = lbcurr(zrcol);
    xzrcol(fubindx) = tub(fubindx);
    xzrcol(flbindx) = tlb(flbindx);
    
    c0curr = c0curr + ccurr(zrcol)' * xzrcol;
    
    nzcol = ~zrcol;
    inzcol = find(nzcol);
    
    if computeLambda 
        % Assign part of Lagrange multipliers according to -ccurr(i), leave others zero
        nn = min(n_orig, ncurr);
        lowermultnonzero = ( ccurr(1:nn,1) >= 0 ) & zrcol(1:nn,1); 
        uppermultnonzero = ( ccurr(1:nn,1) < 0 ) &  zrcol(1:nn,1); 
        lambda.lower(lbindex(lowermultnonzero)) = ccurr(lowermultnonzero); 
        lambda.upper(ubindex(uppermultnonzero)) = -ccurr(uppermultnonzero);
        lbindex = lbindex(inzcol);       %                 
        ubindex = ubindex(inzcol);
        zindex = zindex(inzcol);
        windex = windex(inzcol);
    end
    
    % Update
    Acurr = Acurr(:,nzcol);
    ccurr = ccurr(nzcol);
    lbcurr = lbcurr(nzcol);                        
    ubcurr = ubcurr(nzcol);
    if (verbosity >= 4)
        disp( sprintf('  Deleting %i all-zero columns from the constraint matrix.', ...
            nnz(zrcol)) );
    end
    ncurr = size(ccurr,1);
    colidx = find(MarkVars);
    xpresol(colidx(zrcol)) = xzrcol;
    MarkVars(colidx(zrcol)) = false;   
    
    % Record the recovery info
    restoreSol(2).exist = true;
    restoreSol(2).unslvd_idx = inzcol;  
    restoreSol(2).slvd_idx = izcol;      
    restoreSol(2).slvd_val = xzrcol;
end

% Case: solve singleton rows.
% Here zero inequality rows--which at this stage will have exactly one nonzero
% coefficient corresponding to the slack--will be picked up and considered for
% deletion.
RowSgtonsExist = 0;
rnnzct = sum(spones(Acurr), 2);
singletons = (rnnzct == 1);     
nsgrows = nnz(singletons); 
if (nsgrows >= max(1, .01*size(Acurr,1)))
    RowSgtonsExist = 1;
    nonSgtons = ~singletons; 
    if (verbosity >= 4)
        disp( sprintf('  Solving %i row-singleton variables immediately.', nsgrows) );
        disp('  The row_singleton index:');
        disp( sprintf('   %d', find(singletons)) );
    end
    Atmp  = Acurr(singletons,:);
    Atmp1 = spones(Atmp);
    btmp  = bcurr(singletons);
    if (nsgrows == 1) 
        isolved  = logical(Atmp1);
        unsolved = ~ isolved;
        xsolved  = btmp/Atmp(isolved);
    else
        colnnzct = sum(Atmp1, 1);
        isolved  = logical(spones(colnnzct));
        unsolved = ~ isolved;
        [ii, jj, vv] = find(Atmp); 
        btmp = btmp(ii);  
        xsolved = btmp./vv;
        if ( any(colnnzct > 1) )  % if there is repeating columns in all of the singleton rows
            repeat = diff([0; jj]) == 0; 
            for i = 1: ( length(xsolved) - 1 )
                if repeat(i+1) && ( xsolved(i+1) ~= xsolved(i) )
                    if (verbosity >= 4)
                        disp( sprintf(['  Exiting due to infeasibility:' ...
                                ' Singleton variables in equality constraints are not consistent.']) );
                    end
                    lambda.ineqlin = [];
                    lambda.eqlin = [];
                    lambda.upper = [];
                    lambda.lower = [];
                    lamindx.Aindex = [];
                    PresolExit = -1; % Will be remapped to -2 in linprog.m
                    return
                end
            end % for i = 1: ( length(xsolved) - 1 )
            xsolved(repeat) = [];
        end % if any(colnnzct > 1)
    end % if (nsgrows == 1)
    
    % Check that singleton variables are within bounds
    if ( (any(xsolved < lbcurr(isolved))) || ...
            (any(xsolved > ubcurr(isolved))) )
        if (verbosity >= 4)
            disp( sprintf(['  Exiting due to infeasibility:' ...
                    ' %i singleton variables in the equality constraints are not within bounds.'], ...
                sum((xsolved<lbcurr(isolved))|(xsolved>ubcurr(isolved)))) );
        end
        lambda.ineqlin = [];
        lambda.eqlin = [];
        lambda.upper = [];
        lambda.lower = [];
        lamindx.Aindex = [];
        PresolExit = -1; % Will be remapped to -2 in linprog.m
        return
    end
    
    if computeLambda
        % Compute which Lagrange multipliers need to be computed.
        %   sgrows: what eqlin lambdas we are solving for (row in A)
        %   sgcols: what column in A that singleton is in
        sgAindex = Aindex(singletons);
        sgrows = sgAindex( (sgAindex >= Aeqstart & sgAindex <= Aeqend) ) ...
            - numRowsAineq;
        
        % Only want to extract out indices for the original variables, not slacks.
        % Must also subtract out the removed variables from zero columns and fixed variables.
        [iii, tmp] = find(Atmp1');
        sgcols = lbindex( iii( iii <= (n_orig - nnz(zrcol) - nnz(fixed)) ) );
        if length(sgrows) ~= length(sgcols)
            error(message('optim:simplexpresolve:SizeMismatch'));
        end
        Aindex = Aindex(nonSgtons);
        yindex = yindex(nonSgtons);
        lbindex = lbindex(unsolved);
        ubindex = ubindex(unsolved);
        zindex = zindex(unsolved);
        windex = windex(unsolved);
        
        restoreSol(5).lorow = sgrows;
        restoreSol(5).upcol = sgcols;
    end
    
    % Update
    AindexEq = AindexEq(nonSgtons);
    xpresol(isolved)  = xsolved; 
    c0curr = c0curr + ccurr(isolved)' * xsolved;
    ccurr = ccurr(unsolved);
    bcurr = bcurr(nonSgtons,1) - Acurr(nonSgtons,isolved)*xsolved;
    Acurr = Acurr(nonSgtons, unsolved);
    lbcurr = lbcurr(unsolved);
    ubcurr = ubcurr(unsolved);
    
    mcurr = mcurr - sum(singletons); 
    ncurr = ncurr - nnz(isolved);
    colidx = find(MarkVars);
    MarkVars(colidx(isolved)) = false; 
    if restoreBasis
        rowidx = find(MarkConstr);
        MarkConstr(rowidx(singletons)) = false;
        restoreSol(5).col = colidx(isolved);
        restoreSol(5).row = rowidx(singletons);
    end
    
    % Record the recovery info
    restoreSol(5).exist = true;
    restoreSol(5).unslvd_idx = unsolved;
    restoreSol(5).slvd_idx = isolved;
    restoreSol(5).slvd_val = xsolved; 
end  % if (nsgrows >= max(1, .01*size(Acurr,1)))

% Case: remove forcing constraints
%       compute implied upper bound Rub and implied lower bound Rlb for
%       each row with nonzero elements corresponds to variables with finite lb and ub
if (~isempty(Acurr))
    ubcurr = ubcurr(:);
    lbcurr = lbcurr(:);
    
    Rub = max(Acurr, 0)*ubcurr + min(Acurr, 0)*lbcurr;
    Rlb = max(Acurr, 0)*lbcurr + min(Acurr, 0)*ubcurr;
    
    if (~restoreBasis) 
        % An infeasible constraint detected by incompatibility between bcurr, Rlb & Rub.
        if any( Rlb > bcurr) || any( Rub < bcurr)
            lambda.ineqlin = [];
            lambda.eqlin = [];
            lambda.upper = [];
            lambda.lower = [];
            lamindx.Aindex = [];
            PresolExit = -1; % Will be remapped to -2 in linprog.m
            if verbosity >= 4
                disp('  Exiting due to infeasibility: Detected in forcing constraint case.');
            end
            return
        end
        
        % When Rlb == Rub (== bcurr) hold for some constraints,
        % detect the case that two forcing constraints force one variable to
        % be fixed at different upper and lower bounds at the same time.
        % Note: Rlb <= bcurr <= Rub.
        fixforconstr = (Rub == Rlb);
        if any(fixforconstr)
            conflictfix = (ubcurr ~= lbcurr) & any(Acurr(fixforconstr, :), 1)';
            if any(conflictfix)
                lambda.ineqlin = [];
                lambda.eqlin = [];
                lambda.upper = [];
                lambda.lower = [];
                lamindx.Aindex = [];
                PresolExit = -1; % Will be remapped to -2 in linprog.m
                if verbosity >= 4
                    disp('  Exiting due to infeasibility: Detected in forcing constraint case.');
                end
                return;
            end
        end
        
        % A forcing constraint leads to the fixing of all variables
        % corresponding to nonzero elements in each forcing constraint row.
        uforconstr = ( Rub == bcurr );
        lforconstr = ( Rlb == bcurr );
        ForConstrExist = any( uforconstr | lforconstr );
        
        % Initialize the variables in case that uforconstr/lforconstr may be empty
        uPosAcurr = zeros(1, ncurr);
        uNegAcurr = zeros(1, ncurr);
        lPosAcurr = zeros(1, ncurr);
        lNegAcurr = zeros(1, ncurr);
        
        % Update together for the above two cases according to each column(variable)
        if ForConstrExist
            
            if any(uforconstr)
                uPosAcurr = any( (Acurr(uforconstr, :) > 0), 1 );
                uNegAcurr = any( (Acurr(uforconstr, :) < 0), 1 );
            end
            if any(lforconstr)
                lPosAcurr = any( (Acurr(lforconstr, :) > 0), 1 );
                lNegAcurr = any( (Acurr(lforconstr, :) < 0), 1 );
            end
            
            % Detect the case that two forcing constraints force one variable to
            % be fixed at different upper and lower bounds at the same time.
            ubfixcands = (uPosAcurr | lNegAcurr )';
            lbfixcands = (uNegAcurr | lPosAcurr )';
            confix = (lbcurr ~= ubcurr) & ubfixcands  & lbfixcands;
            if any(confix)
                lambda.ineqlin = [];
                lambda.eqlin = [];
                lambda.upper = [];
                lambda.lower = [];
                lamindx.Aindex = [];
                PresolExit = -1; % Will be remapped to -2 in linprog.m
                if verbosity >= 4
                    disp('  Exiting due to infeasibility: Detected in forcing constraint case.');
                end
                return;
            end
            
            varcurr = find(MarkVars);
            fub = ubcurr(ubfixcands);
            flb = lbcurr(lbfixcands);
            xpresol(varcurr(ubfixcands)) =  fub;
            xpresol(varcurr(lbfixcands)) =  flb;
            ffixed = ubfixcands | lbfixcands;
            MarkVars( varcurr(ffixed) ) = false;
            
            if restoreBasis
                constrcurr = find(MarkConstr);
                MarkConstr(constrcurr(uforconstr | lforconstr)) = false;
                restoreSol(4).row = constrcurr(uforconstr | lforconstr);
                restoreSol(4).col = varcurr(ffixed);
            end
            
            rowforce = (Rub == bcurr) | (Rlb == bcurr);
            % Keep the unremoved (See presolve Case: remove zero rows) zero rows
            % that are marked as forcing constraints.
            % Check if the relevant Acurr section has a zero row.
            % Relevant Acurr section is the part of Acurr that corresponds to the
            % forcing costraints, and variables that are not yet solved for in presolve, i.e.,
            % find(rowforce), 1:nnz(MarkVars(1:n_orig)).
            indxFRr = find(rowforce)';
            indxFRc = 1:nnz(MarkVars(1:n_orig));
            if ~isempty(Acurr(indxFRr,indxFRc))
                zeroForceRow = (sum(spones(Acurr(indxFRr,indxFRc)),2) == 0);
                rowforce(indxFRr(zeroForceRow)) = false;
            end
            
            rownf = ~rowforce;
            nfix = ~ffixed;
            
            if verbosity >= 4
                disp( sprintf(['  Fixing %d variables at upper bounds and %d variables at lower bounds: ' ...
                    'Due to %d forcing constraints.'], nnz(ubfixcands), nnz(lbfixcands), nnz(rowforce) ));
            end
            
            if computeLambda
                restoreSol(4).lorow = lbindex( lbfixcands & (lbindex' <= n_orig) );
                restoreSol(4).upcol = ubindex( ubfixcands & (ubindex' <= n_orig) );
                restoreSol(4).Adel = Aindex( rowforce );
                restoreSol(4).bdel = lforconstr(  lforconstr|uforconstr );
                
                % Update index
                Aindex = Aindex(rownf);
                yindex = yindex(rownf);
                lbindex = lbindex(nfix);
                zindex = zindex(nfix);
                ubindex = ubindex(nfix);
                windex = windex(nfix);
            end % end of "if computeLambda"
            
            % Update
            AindexEq = AindexEq(rownf);
            c0curr = c0curr  + [ccurr(ubfixcands)' ccurr(lbfixcands)'] * [fub; flb];
            ccurr  = ccurr(nfix);
            bcurr  = bcurr(rownf, 1) - [Acurr(rownf, ubfixcands) Acurr(rownf, lbfixcands)] * [fub;flb];
            Acurr  = Acurr(rownf, nfix);
            lbcurr = lbcurr(nfix);
            ubcurr = ubcurr(nfix);
            [mcurr, ncurr]  = size(Acurr);
            
            % Update Rub and Rlb for their usage in dominated constraint
            % procedure to detect implied free variables.
            Rub = Rub(rownf);
            Rlb = Rlb(rownf);
            
            % Record the recovery info
            restoreSol(4).exist      = true;
            restoreSol(4).unslvd_idx = nfix;
            restoreSol(4).slvd_idx   = ffixed;
            xtmp = zeros(length(ffixed),1);
            xtmp(ubfixcands, 1)  = fub;
            xtmp(lbfixcands, 1)  = flb;
            restoreSol(4).slvd_val   = xtmp(ffixed);
        end % end of "if ForConstrExist"
    end % if ~restoreBasis
end % end if ~isempty(Acurr) forcing constraint case

% Case: dominated constraints (together with forcing constraints)
%       compute implied bounds on variable xj
%       combine these bounds with the original bounds to detect more free column singletons
%       only compute the implied bounds for each aij != 0 if either Rub or Rlb are finite 

% Limit on detecting dominated constraints because of the computational cost.
if (~isempty(Acurr)) && ( nnz(Acurr)/(mcurr*ncurr) < 0.1 )
    nnzAcurr = nnz(Acurr);
    ubImplied = inf(nnzAcurr, 1);
    lbImplied = -ubImplied;
    
    [iarray, jarray, varray] = find(Acurr); 
    iarray= iarray(:); 
    jarray = jarray(:);
    varray = varray(:); 
    lbij = lbcurr(jarray);
    ubij = ubcurr(jarray); 
    
    Rlbij = Rlb(iarray);
    Rubij = Rub(iarray);
    bij = bcurr(iarray);
    
    indposAij =  find( varray > 0 ); % the index of arrays: iarray, jarray and varray 
    indnegAij =  find( varray < 0 );                                                 
    
    if ~isempty( indposAij )
        % The index of the subarray Rubij(indposAij)
        pindfiniteRubij = isfinite( Rubij(indposAij) );             
        pindfiniteRlbij = isfinite( Rlbij(indposAij) );             
        if any( pindfiniteRubij )    
            % Transform back to the index of iarray.
            pfiniteRubij = indposAij(pindfiniteRubij);    
            lbImplied(pfiniteRubij) = ubij(pfiniteRubij) + ...
                ( bij(pfiniteRubij) - Rubij(pfiniteRubij) ) ./ varray(pfiniteRubij);
        end  
        if any( pindfiniteRlbij ) 
            pfiniteRlbij = indposAij(pindfiniteRlbij); 
            ubImplied(pfiniteRlbij) = lbij(pfiniteRlbij) + ...
                ( bij(pfiniteRlbij) - Rlbij(pfiniteRlbij) ) ./ varray(pfiniteRlbij);
        end
    end
    if ~isempty( indnegAij )
        nindfiniteRubij =  isfinite( Rubij(indnegAij) );        
        nindfiniteRlbij =  isfinite( Rlbij(indnegAij) );       
        if any( nindfiniteRlbij )
            nfiniteRlbij = indnegAij(nindfiniteRlbij); 
            lbImplied(nfiniteRlbij) = ubij(nfiniteRlbij) + ... 
                ( bij(nfiniteRlbij) - Rlbij(nfiniteRlbij) ) ./ varray(nfiniteRlbij);
        end
        if any( nindfiniteRubij )
            nfiniteRubij = indnegAij(nindfiniteRubij);
            ubImplied(nfiniteRubij) = lbij(nfiniteRubij) + ...
                ( bij(nfiniteRubij) - Rubij(nfiniteRubij) ) ./ varray(nfiniteRubij);
        end                                                       
    end
    
    % Find locations of unique indices in jarray (which is sorted)
    groups = find(diff([0; jarray(:)]) ~= 0);
    ubImp = inf(ncurr, 1);
    lbImp = -ubImp;
    ubImp(jarray(groups)) = ubImplied(groups);
    lbImp(jarray(groups)) = lbImplied(groups);
    ngroups = length(groups);
    
    % If there exists repeating implied bounds on the same variable
    if ngroups < ncurr 
        for k = 1 : (ngroups-1)
            ubImp(jarray(groups(k))) = min( ubImplied(groups(k) : (groups(k+1)-1)) );
            lbImp(jarray(groups(k))) = max( lbImplied(groups(k) : (groups(k+1)-1)) );
        end
        k = groups(ngroups);  % k is the beginning index of the last group(column)
        if k < nnzAcurr
            ubImp( jarray(k) ) = min( ubImplied(k: nnzAcurr) );   
            lbImp( jarray(k) ) = max( lbImplied(k: nnzAcurr) );
        end
    end
    freeImp = (lbcurr <= lbImp) & (ubImp <= ubcurr);
    %------------------------------------------------------------
    
    % case: Remove all free and implied free column singletons
    % colcurr has 3 values: 0, not a free/implied-free column singleton; 
    %                       1, free/implied-free column singleton with duplicate column; 
    %                      -1, to be deleted free/implied-free column singleton
    colcurr = zeros(length(lbcurr), 1); 
    % Find column singletons: columns that have exactly one non-zero
    % coefficient in Acurr, and this non-zero coefficient corresponds 
    % to an equality. Also, to be a column singleton, the associated
    % variable has to be free (or implied free).
    sumcolnz   = sum( spones(Acurr),1 );             % # of nonzero coeff's per col in Acurr
    sumcolnzEq = sum( spones(Acurr(AindexEq,:)),1 ); % # of nonzero coeff's per col in equalities
    colcurr = (sumcolnzEq == 1)' & (sumcolnz == 1)' & ( (isinf(ubcurr) & isinf(lbcurr)) | freeImp );    
        
    indcolSgtons = find( colcurr );                  
    ColSgltonsExist = ~isempty( indcolSgtons ); 
    % To avoid the case that the deletion of singleton columns may lead to
    % delete all rows, add the size condition.
    if  ColSgltonsExist && ( length(bcurr) > nnz(colcurr) )
        subAcurr= Acurr(:, indcolSgtons); 
        sumrsubAcurr = sum( spones(subAcurr), 2 );
        rowdel =  (sumrsubAcurr == 1);
        if nnz(rowdel) < length(rowdel)
            rowleft = ~ rowdel; 
        else 
            rowleft = [];             
        end
        % Note: the case that sumrsubAcurr >= 2, a simple case of multiple columns.
        
        ssubAcurr = subAcurr(rowdel, :);         
        [csrowdel, coldel, valdel] = find(ssubAcurr);                            
        coldel = indcolSgtons(coldel);           
        colcurr = double (colcurr);
        colcurr(coldel) = -1;        % -1 represents that the column is deleted. 
        colleft =  ( colcurr >= 0 );             
        
        % Update the objective function when necessary ------------
        if any(rowdel) && (~isempty(coldel))    
            if (verbosity >= 4)
                disp( sprintf('  Deleting  %i free or implied free singleton columns.', ...
                    nnz(coldel)) );
            end
            % Getting the order of the deleted rows wrt deleted columns 
            % For each column singleton we remove the corresponding column and row from the problem.
            % In postsolve these deleted rows are used to calculate the values of the variables 
            % that are eliminated. The columns and rows that are eliminated are sorted in ascending order, 
            % however the calculations for the eliminated columns should be done considering their matching row.
            % For example, assume that we have four column singletons, namely variables 5,12,13,15 
            % and the corresponding rows are 215,25,218,212.
            %
            %            5 ... 12 13 ... 15
            %      ...      
            %       25          x             
            %      ...
            %      212                    x
            %      ... 
            %      215   x
            %      ...
            %      218             x
            %
            % However, the vector rowdel does not contain this ordering information, instead its usage
            % results in a deleted-row ordering, 25,212,215,218, or equivalently (in terms of deleted rows)
            % 1,2,3,4. The correct ordering should be 215,25,218,212 (or equivalently 3,1,4,2).
            %
            % The correct ordered row indices can be obtained by
            % find(subAcurr) - (0:mcurr:mcurr*(length(indcolSgtons)-1))'
            % The same ordering information with an offset is also available from vector csrowdel above.
            % To transform the ordering 215,25,218,212 into 3,1,4,2 we use the sort function twice. The second
            % output argument of the first sort gives the original location indices corresponding to the 
            % sorted vector, namely 2,4,1,3. The second output argument of the second sort called with this vector 
            % results in the ordering that we are after, namely 3,1,4,2. 
            
            [dummysort,ordrowdelindx] = sort(csrowdel);
            [dummysort,ordrowdelindx] = sort(ordrowdelindx);
            
            if computeLambda
                % Record the information on the deleted rows and columns
                restoreSol(3).lorow = Aindex(rowdel) - Aineqend;  % indices w.r.t. Aeq
                % Ordering the rows according to corresponding singleton  
                % columns. Here restoreSol(3).lorow is used in postsolve to calculate 
                % the lambda's corresponding to the deleted rows. These rows are deleted 
                % because they contain a column singleton.
                restoreSol(3).lorow = restoreSol(3).lorow(ordrowdelindx);
                
                restoreSol(3).upcol = lbindex(coldel);
                Aindex = Aindex(rowleft);
                yindex = yindex(rowleft);
                % Need to make sure the deleted row are in the equality constraint 
                % lambda.eqlin(rowdel) = ccurr(coldel)./valdel.
                lbindex = lbindex(colleft);
                zindex = zindex(colleft);
                ubindex = ubindex(colleft);
                windex = windex(colleft);
            end
            
            % Update           
            AindexEq = AindexEq(rowleft);
            bdel = bcurr(rowdel);
            Adel = Acurr(rowdel, colleft);
            % Ordering the rows according to corresponding singleton
            % columns
            bdel = bdel(ordrowdelindx);
            Adel = Adel(ordrowdelindx,:);
            
            c0curr = c0curr + sum( ccurr(coldel) ./ valdel .* bdel ) ;    
            ccurr  = ccurr(colleft) - Adel' * (ccurr(coldel)./valdel);
            
            Acurr  = Acurr(rowleft, colleft);  
            bcurr  = bcurr(rowleft);
            lbcurr = lbcurr(colleft);
            ubcurr = ubcurr(colleft);
            [mcurr, ncurr] = size(Acurr);
            
            colidx = find(MarkVars);
            MarkVars(colidx(coldel)) = false;
            if restoreBasis
                rowidx = find(MarkConstr);
                MarkConstr(rowidx(rowdel)) = false;
                restoreSol(3).row = rowidx(rowdel);
                restoreSol(3).col = colidx(coldel);
            end
            
            % Record the recovery info
            restoreSol(3).exist = true;
            restoreSol(3).unslvd_idx = colleft;
            restoreSol(3).slvd_idx = coldel;
            restoreSol(3).slvd_val = valdel;
            restoreSol(3).Adel = Adel;
            restoreSol(3).bdel = bdel;
        end % end of "if any(rowdel) && (~isempty(coldel))" 
        
    end % end of "if  ColSgltonsExist && ( length(bcurr) > nnz(colcurr) )"
    
    % It is possible that new zero columns are generated in this case.
    % Add the check for zero columns. 
    ZrcolsExist2 = 0;
    if isempty(Acurr)
        zrcol = 0;
    else
        zrcol = (max(abs(Acurr)) == 0)';  %max behaves differently if A only has one row:
                                          %so zero cols won't be deleted in this case.
    end
    
    if ( any(zrcol == 1) )
        ZrcolsExist2 = 1;
        izcol = find(zrcol);
        if ( any( (ccurr(zrcol) < 0) & isinf(ubcurr(zrcol)) ) ...
                ||  any( (ccurr(zrcol) > 0) & isinf(lbcurr(zrcol)) ) )
            if (verbosity >= 4)
                disp( sprintf(['  Exiting due to unboundedness: The problem ',... 
                               'is unbounded detected by generated zero columns.']) );
            end
            lambda.ineqlin = [];
            lambda.eqlin = [];
            lambda.upper = [];
            lambda.lower = [];
            lamindx.Aindex = [];
            PresolExit = -2; % Will be remapped to -3 in linprog.m
            return;
        end
        xzrcol = zeros(nnz(zrcol), 1);
        fubindx = ccurr(zrcol) < 0;
        flbindx = ccurr(zrcol) > 0;
        tub = ubcurr(zrcol); 
        tlb = lbcurr(zrcol);
        xzrcol(fubindx) =  tub(fubindx);
        xzrcol(flbindx) =  tlb(flbindx);
        
        c0curr = c0curr + ccurr(zrcol)' * xzrcol;
        
        nzcol = ~zrcol;
        inzcol = find(nzcol);
        
        if computeLambda
            % Assign part of Lagrange multipliers according to -ccurr(i), leave others zero 
            lowermultnnz = ( ccurr(1:n_orig,1) >= 0 ) & zrcol(1:n_orig,1); 
            uppermultnnz = ( ccurr(1:n_orig,1) < 0 ) &  zrcol(1:n_orig,1);  
            lambda.lower(lbindex(lowermultnnz)) = ccurr(lowermultnnz); 
            lambda.upper(ubindex(uppermultnnz)) = -ccurr(uppermultnnz);
            lbindex = lbindex(inzcol);                        
            ubindex = ubindex(inzcol);
            zindex = zindex(inzcol);
            windex = windex(inzcol);
        end
        
        % Update 
        Acurr = Acurr(:,nzcol);
        ccurr = ccurr(nzcol);
        lbcurr = lbcurr(nzcol);                        
        ubcurr = ubcurr(nzcol);
        if (verbosity >= 4)
            disp( sprintf('  Deleting %d generated all-zero columns from constraint matrix.', ...
                nnz(zrcol)) );
        end
        [mcurr, ncurr] = size(Acurr);
        colidx = find(MarkVars);
        xpresol(colidx(zrcol)) = xzrcol;
        MarkVars(colidx(zrcol)) = false;   
        
        % Record the recovery info
        restoreSol(7).exist = true;
        restoreSol(7).unslvd_idx = inzcol;  
        restoreSol(7).slvd_idx = izcol;      
        restoreSol(7).slvd_val = xzrcol;
    end
    
end % end of if ~isempty(Acurr) before the dominated constraint procedure

% Case: delete zero rows
ZrowsExist = 0;
rnnzct = sum(spones(Acurr'), 1);     
if ~isempty(Acurr) && (any(rnnzct == 0))
    zrows = (rnnzct == 0);
    ZrowsExist = any(zrows);
    if (any(bcurr(zrows) ~= 0))
        if (verbosity >= 4)
            disp( sprintf(['  Exiting due to infeasibility: An all-zero row in the constraint' ...
                    ' matrix does not have a zero in corresponding right-hand-side entry.']) );
        end
        lambda.ineqlin = [];
        lambda.eqlin = [];
        lambda.upper = [];
        lambda.lower = [];
        lamindx.Aindex = [];
        PresolExit = -1; % Will be remapped to -2 in linprog.m
        return
    else   
        nzrows = ~zrows;
        Acurr  = Acurr(nzrows,:);
        bcurr  = bcurr(nzrows); 
        if (verbosity >= 4)
            disp( sprintf('  Deleting %i all-zero rows from the constraint matrix.', nnz(zrows)) );
        end
        
        if computeLambda
            % Lagrange multipliers for zero rows are zero
            Aindex = Aindex(nzrows);
            yindex = yindex(nzrows);
        end
        AindexEq = AindexEq(nzrows);            
        
        if restoreBasis
            rowidx = find(MarkConstr);
            MarkConstr(rowidx(zrows)) = false;
        end
        mcurr = size(bcurr,1);
    end
end

% Case: row and column scaling 
% Assumption: no zero row or zero column exists.
badscl = 1.e-4;
ColScaled = 0;
absnzs = abs( nonzeros(Acurr) );
thescl = min(absnzs) / max(absnzs);
Ubounds_exist = full(any(ubcurr ~= Inf)); 
Lbounds_exist = full(any(lbcurr ~= -Inf));

if (thescl < badscl)
    if (verbosity >= 4)
        disp(['  Scaling problem by square roots of infinity norms of rows and' ...
                ' columns of constraint matrix.']);
    end
    
    % ----- Scaling vectors ------
    absA = abs(Acurr);
    colscl = full( sqrt(max(absA, [], 1)') );
    rowscl = full( sqrt(max(absA,[], 2)) );
        
    % ----- Column scaling -----
    if (Ubounds_exist)
        ubcurr = ubcurr .* colscl;
    end
    
    if (Lbounds_exist)
        lbcurr = lbcurr .* colscl; % necessary as lb doesn't have to be zero
    end
    
    colscl = reciprocal(colscl); 
    Acurr = Acurr * spdiags(colscl,0,ncurr,ncurr);
    
    ccurr = ccurr .* colscl;
    ColScaled = 1;
    
    % ----- Row scaling -----
    rowscl = reciprocal(rowscl);
    Acurr = spdiags(rowscl,0,mcurr,mcurr) * Acurr;
    bcurr = bcurr .* rowscl;
    bnrm = norm(bcurr);
    q = 1;
    if (bnrm > eps) 
        q = median([1 norm(ccurr)/bnrm 1.e+8]);
        if (q > 10)
            Acurr = q * Acurr;
            bcurr = q * bcurr;
        else
            q = 1;
        end
    end
     
    % Record the recovery info
    restoreSol(1).exist = true;
    restoreSol(1).slvd_val = colscl;
    restoreSol(1).upcol = q;
    restoreSol(1).lorow = rowscl;
end % end of "if (thescl < badscl)"

% Check whether all constraints and variables have been deleted from previous presolve.
if isempty(Acurr) 
    PresolExit = 0;
end

if isempty(Acurr) && ~isempty(lbcurr)
    % Only the lb anf ub boundary conditions need to be satisfied:
    % fix variables at their bound according to their objective coefficients.
    
    nccurr = ccurr <= 0;
    pccurr = ccurr > 0;
    xu = ubcurr(nccurr);
    xl = lbcurr(pccurr);
    xfix(pccurr,1)= xl;
    xfix(nccurr,1)= xu;
    
    if nnz(isinf(xfix) > 0 )
        lambda.ineqlin = [];
        lambda.eqlin = [];
        lambda.lower = [];
        lambda.upper = [];
        lamindx.Aindex = [];
        PresolExit = -2; % Will be remapped to -3 in linprog.m
        if verbosity >= 4
            disp( sprintf(['  Exiting due to unboundedness: The problem is ',...
                  'unbounded due to fixing variables at infinite bounds.']) );
        end
        return;
    end
    
    if computeLambda
        lambda.lower(lbindex(pccurr)) = ccurr(pccurr);
        lambda.upper(ubindex(nccurr)) = -ccurr(nccurr);
        lbindex = [];
        ubindex = [];
        zindex = [];
        windex = [];
    end
    
    restoreSol(8).exist = true;
    restoreSol(8).unslvd_idx = [];
    restoreSol(8).slvd_idx = 1:ncurr;   
    restoreSol(8).slvd_val = xfix;
    colidx = find(MarkVars);
    MarkVars(colidx) = false;
    restoreSol(8).col = colidx; 
    c0curr = c0curr + ccurr'*xfix;
    
    if verbosity >=4 
        disp('  The constraint matrix becomes empty after preprocessing.'); 
    end
end % end of "if isempty(Acurr) && ~isempty(lbcurr)"

ndel = n_orig - nnz(MarkVars(1:n_orig));
% Collect all the index information into a structure lamindx.
if computeLambda
    lamindx.Aindex = Aindex;
    lamindx.yindex = yindex;
    lamindx.lbindex = lbindex;
    lamindx.ubindex = ubindex;
    lamindx.zindex = zindex;
    lamindx.windex = windex;
    
    lamindx.lbend = lbend;
    lamindx.ubend = ubend;
    lamindx.Aineqstart = Aineqstart;
    lamindx.Aineqend = Aineqend;
    lamindx.Aeqstart = Aeqstart;
    lamindx.Aeqend = Aeqend;
end

% End of the SIMPLEXPRESOLVE procedure

%--------------------------------------------------------------------------
%------------------------SIMPLEXPRESOLVE Subfunction-----------------------
%--------------------------------------------------------------------------

function Y = reciprocal(X)
%RECIPROCAL  Invert the nonzero entries of a matrix elementwise.
%   Y = RECIPROCAL(X) has the same sparsity pattern as X
%	 (except possibly for underflow).

if (issparse(X))
    [m, n]  = size(X);
    [i,j,Y] = find(X);
    Y = sparse(i,j,1./Y,m,n);
else
    Y = 1./X;
end

Y = min(Y,1e8); 

%--------------------------------------------------------------------------
