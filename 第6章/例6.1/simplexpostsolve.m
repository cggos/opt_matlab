function [xsol, fval, lambda, exitflagps, basicVarIdx, nonbasicVarIdx, ...
    delrows] = simplexpostsolve(x, fv, c0, c_orig, restoreSol, ...
    diagnostic_level, lambda, lamindx, dualvars, Aineq_orig, Aeq_orig, ...
    computeLambda, caseflag, MarkVars, rbasicVarIdx, rnonbasicVarIdx, ...
    restoreBasis, MarkConstr, delrows1)
%SIMPLEXPOSTSOLVE Restore the solution and the Lagrange multipliers lambda.

%   Copyright 1990-2014 The MathWorks, Inc.

% Initialization
exitflagps = 1;
xsol = x;
n = length(xsol);
n_orig = length(c_orig);
xtmp = zeros(n_orig + size(Aineq_orig, 1), 1);

% Restore the solution to the original system
% Case: fix variables at bounds according to the objective coefficient 
if (restoreSol(8).exist)
    xtmp(restoreSol(8).unslvd_idx,1) = xsol;
    fixed = restoreSol(8).slvd_idx;
    xtmp(fixed, 1) = restoreSol(8).slvd_val;
    n = length(fixed) + length(xsol);
    xsol = xtmp(1:n);
    if (diagnostic_level >= 4)
        disp( sprintf('  Restoring %d fixed variables with solution.', ...
            nnz(fixed)) );
    end
end

% Case: column scaling
if (restoreSol(1).exist)
    xsol = xsol .* restoreSol(1).slvd_val(1:n);
    if (diagnostic_level >= 4)
        disp('  Scaling solution back to original problem.');
    end 
end

% Case: generated zero column case 
if (restoreSol(7).exist)
    xtmp(restoreSol(7).unslvd_idx,1) = xsol;
    xzrcol = restoreSol(7).slvd_val;
    xtmp(restoreSol(7).slvd_idx,1) = xzrcol;
    n = length(xsol) + length(xzrcol);
    xsol = xtmp(1:n);
    if (diagnostic_level >= 4)
        disp( sprintf('  Restoring %d generated zero variables with solution.', ...
            length(xzrcol)) );
    end
end

% Case: free and implied free column singleton
% free column singleton: A(i, :), b(i)  where a(i,j)~= 0 and column j is a column singleton
% data needed: Adel = A(rowdel, colleft), bdel, valdel for deleted rows
%                     coldel, the index of the deleted columns
%                     colleft, for the remaining columns
%                     to compute the value of the deleted singleton columns
if (restoreSol(3).exist)
    xtmp(restoreSol(3).unslvd_idx,1) = xsol;
    coldel = restoreSol(3).slvd_idx;
    xtmp(coldel,1) = ( restoreSol(3).bdel - restoreSol(3).Adel * xsol) ./  restoreSol(3).slvd_val;
    n = length(restoreSol(3).unslvd_idx);
    xsol = xtmp(1:n);
    if (diagnostic_level >= 4)
        disp( sprintf('  Restoring %d free or implied free column singleton variables with solution.', ...
            length(coldel)) );
    end
end

% Case: forcing constraints
if (restoreSol(4).exist)
    xtmp(restoreSol(4).slvd_idx,1) = restoreSol(4).slvd_val; 
    nfix = restoreSol(4).unslvd_idx;
    xtmp(nfix,1) = xsol;
    n = length(nfix);
    xsol = xtmp(1:n);
    if (diagnostic_level >= 4)
        disp( sprintf('  Restoring %d fixed variables by forcing constraints.', ...
            nnz(~nfix)) );
    end
end

% Case: row singletons
if (restoreSol(5).exist)
    xtmp(restoreSol(5).unslvd_idx,1) = xsol;
    xsolved = restoreSol(5).slvd_val;
    xtmp(restoreSol(5).slvd_idx,1) = xsolved;
    n = length(restoreSol(5).slvd_idx);
    xsol = xtmp(1:n); 
    if (diagnostic_level >= 4)
        disp( sprintf('  Restoring %d row singleton variables with solution.', ...
            length(xsolved)) );
    end
end

% Case: zero column 
if (restoreSol(2).exist)
    xtmp(restoreSol(2).unslvd_idx,1) = xsol;
    xzrcol = restoreSol(2).slvd_val;
    xtmp(restoreSol(2).slvd_idx,1) = xzrcol;
    n = length(xsol) + length(xzrcol);
    xsol = xtmp(1:n);
    if (diagnostic_level >= 4)
        disp( sprintf('  Restoring %d obvious zero variables with solution.', ...
            length(xzrcol)) );
    end
end

% Case: dependent rows deleted
%       need to check the consistency for feasibility
if (restoreSol(9).exist)
    Adel = restoreSol(9).Adel;
    bdel = restoreSol(9).bdel;
    if (norm(Adel*xsol-bdel)/max(1,norm(bdel)) > 1.0e-8)     
        exitflagps = -1; % Will be remapped to -2 in linprog.m
        if (diagnostic_level >= 4)
            disp(sprintf(['  The primal problem is infeasible: \n' ...
                    '    Some equality constraints are dependent but not consistent.']));
        end
    else
        if (diagnostic_level >= 4)
            disp( sprintf('  Solution satisfies %d dependent constraints.', ...
                size(Adel,1)) );
        end
    end
end

% Case: fixed variables by equal upper and lower bounds
if (restoreSol(6).exist)
    xtmp(restoreSol(6).unslvd_idx,1) = xsol;
    fixed = restoreSol(6).slvd_idx;
    xtmp(fixed, 1) = restoreSol(6).slvd_val;
    n = length(fixed);
    xsol = xtmp(1:n);
    if (diagnostic_level >= 4)
        disp( sprintf('  Restoring %d fixed variables (due to equal bounds) with solution.', ...
            nnz(fixed)) );
    end
end

% Case: slack variables added
if  n > n_orig   
    % Remove slack solution variables
    xsol = xsol(1:n_orig);
    if (diagnostic_level >= 4)
        disp( sprintf('  Eliminating %d slack variables from solution.', ...
            n-n_orig) );
    end
end

% Compute the objective value: two ways to the same answer 
fval0 = fv + c0; 
fval  = full(c_orig' * xsol);

if  diagnostic_level >=4
    % for debug purpose
    if ( (abs(fval0) > 1.0e-12)  && (abs(fval - fval0)/abs(fval0) > 1.e-6) ) || ... % relative error
            ( (abs(fval0) < 1.0e-12) && (abs(fval - fval0) > 1.e-6) )
        if (diagnostic_level >= 4)
            disp( sprintf('  fval = %f, \t fval0 = %f\n', fval, fval0) ); 
        end
        error(message('optim:simplexpostsolve:RelErrors'));
    end
end

%--------------------------------------
% restoreSol(9).case = 'RowDependent';
% restoreSol(8).case = 'FixedVarsbyObj';
% restoreSol(7).case = 'ZeroColsGenerated';
% restoreSol(6).case = 'FixedVars';
% restoreSol(5).case = 'RowSgtons';
% restoreSol(4).case = 'ForConstr';
% restoreSol(3).case = 'ColSgtons';
% restoreSol(2).case = 'ZeroCols';
% restoreSol(1).case = 'ColScal';
%--------------------------------------

% Restore Lambda if required.
if computeLambda
    % Calculate Lagrange multipliers in two cases
    if caseflag == 1
        % Case 1. presolve eliminated all the constraints, restore lambda
        %         from the presolve directly.
        if restoreSol(3).exist % Column Singletons
            % csgrows is the index of the deleted row according to column singletons
            % csgcols is the index of the column singletons
            csgrows = restoreSol(3).lorow;
            csgcols = restoreSol(3).upcol;
            lambda.eqlin(csgrows) =  -(c_orig(csgcols))./ Aeq_orig(csgrows,csgcols);
        else
            csgrows = []; 
            csgcols = [];
        end
        
        if restoreSol(4).exist % Forcing Constraints
            ffixedlower = restoreSol(4).lorow;
            ffixedupper = restoreSol(4).upcol;
            Rlbeqb = restoreSol(4).bdel; % Rubeqtag = ~Rlbeqtag;
            fixforce = [ffixedlower; ffixedupper];
            
            indxforcerows = restoreSol(4).Adel;
            nforcerows = length(indxforcerows);
            frAineqindx = indxforcerows(indxforcerows >= lamindx.Aineqstart ...
                & indxforcerows <= lamindx.Aineqend);
            frAeqindx = indxforcerows(indxforcerows >= lamindx.Aeqstart ...
                & indxforcerows <= lamindx.Aeqend) - size(Aineq_orig, 1);
            
            rc = c_orig(fixforce);
            if ~isempty(csgrows)
                rc = Aeq_orig(csgrows, fixforce)' * lambda.eqlin(csgrows) + rc;
            end
             
             
            A = [Aineq_orig(frAineqindx, fixforce);
                Aeq_orig(frAeqindx, fixforce)];
            
            yforce = zeros(nforcerows, 1);
            for i = 1: nforcerows
                tcol = A(i, :) ~= 0;
                if Rlbeqb(i) % in the case Rlb_i = b_i
                    [ty, tindx] = min( rc(tcol)./ A(i, tcol)' );
                else % in the case Rub_i = b_i  
                    [ty, tindx] = max( rc(tcol)./ A(i, tcol)' );
                end
                yforce(i) = ty;
                rc(tcol) = rc(tcol) - ty * A(i, tcol)';
            end
            
            lambda.ineqlin(frAineqindx) = -yforce(1:length(frAineqindx));
            lambda.eqlin(frAeqindx) = -yforce((length(frAineqindx)+1):end);
            
            lambda.lower(ffixedlower) = rc(1:length(ffixedlower));
            lambda.upper(ffixedupper) = -rc((length(ffixedlower)+1):end); 
        end % if restoreSol(4).exist
        
        if restoreSol(5).exist % Row Singletons
            % sgrows: what eqlin lambdas we are solving for (row in A)
            % sgcols: what column in A that singleton is in
            sgrows = restoreSol(5).lorow;
            sgcols = restoreSol(5).upcol;
            result = -c_orig + lambda.lower - lambda.upper;
            if ~isempty(sgrows)
                lambda.eqlin(sgrows) = result(sgcols) ./ diag(Aeq_orig(sgrows, sgcols)); 
            end
        end
        
        if restoreSol(6).exist % Fixed variables by equal bounds
            fixedlower = restoreSol(6).lorow;
            fixedupper = restoreSol(6).upcol;
            Alambda = Aineq_orig'*lambda.ineqlin + Aeq_orig'*lambda.eqlin;
            lambda.lower(fixedlower) = c_orig(fixedlower) + Alambda(fixedlower);
            lambda.upper(fixedupper) = - c_orig(fixedupper) - Alambda(fixedupper);
        end
        
    else 
        % Case 2. recover lambda from the reduced system and the presolve-----
        % make assignment to lambda.eqlin and lambda.ineqlin.
        % dual variables y, z and w were assigned in function simplexphasetwo.
        lbindex = lamindx.lbindex;
        ubindex = lamindx.ubindex;
        Aindex  = lamindx.Aindex;
        
        lbindexfinal = lbindex(lbindex <= lamindx.lbend); 
        ubindexfinal = ubindex(ubindex <= lamindx.ubend);
        zindexfinal  = find(lamindx.zindex);
        windexfinal  = find(lamindx.windex);
        nnzlbindexfinal = nnz(lbindexfinal);
        nnzubindexfinal = nnz(ubindexfinal);
        
        % If the scaling case happens, 
        % transform the dual variables to the original system
        if restoreSol(1).exist
            colscl = restoreSol(1).slvd_val;
            rowscl = restoreSol(1).lorow;
            alpha  = restoreSol(1).upcol;
            % Use lamindx.yindex to include the case simplexphaseone may delete dependent rows.
            dualvars.y = (dualvars.y) .* rowscl(lamindx.yindex) * alpha; 
            dualvars.z = (dualvars.z) ./ colscl;
            dualvars.w = (dualvars.w) ./ colscl;
        end
        
        y = dualvars.y;
        z = dualvars.z;
        w = dualvars.w;
        
        iAineqindex = Aindex( Aindex >= lamindx.Aineqstart & Aindex <= lamindx.Aineqend );
        % iAineqindex = Aindex( find (Aindex >= lamindx.Aineqstart & Aindex <= lamindx.Aineqend) );
        % Subtract off numRowsAineq even if those were removed since the indices still reflect their existence
        numRowsAineq = size(Aineq_orig,1);
        iAeqindex = Aindex(Aindex >= lamindx.Aeqstart & Aindex <= lamindx.Aeqend) - numRowsAineq;
        
        nnzAineqindex = nnz(iAineqindex);
        nnzAeqindex = nnz(iAeqindex);
        
        lambda.ineqlin(iAineqindex) = full(-y(1:nnzAineqindex));
        yAeq = y(nnzAineqindex+1:end);
        lambda.eqlin(iAeqindex) = full(-yAeq(1:nnzAeqindex));
        
        lambda.upper(ubindexfinal) = full(w(windexfinal(1:nnzubindexfinal)));
        lambda.lower(lbindexfinal) = full(z(zindexfinal(1:nnzlbindexfinal)));
                  
        % --------Restore lambda from the presolve case by case--------
        % Case: column singleton
        if restoreSol(3).exist
            % csgrows is the index of the deleted row according to column singletons
            % csgcols is the index of the column singletons.
            % diag(Aeq_orig(csgrows,csgcols)) contains the column singleton coefficients.
            csgrows = restoreSol(3).lorow;
            csgcols = restoreSol(3).upcol;
            lambda.eqlin(csgrows) =  -(c_orig(csgcols))./ diag(Aeq_orig(csgrows,csgcols));
        else
            csgrows = []; 
            csgcols = [];
        end
              
        % Case: forcing constraints
        if restoreSol(4).exist
            ffixedlower = restoreSol(4).lorow;
            ffixedupper = restoreSol(4).upcol;
            Rlbeqb = restoreSol(4).bdel; 
            fixforce = [ffixedlower(:); ffixedupper(:)];
            
            indxforcerows = restoreSol(4).Adel;
            nforcerows = length(indxforcerows);
            frAineqindx = indxforcerows(indxforcerows >= lamindx.Aineqstart ...
                & indxforcerows <= lamindx.Aineqend);
            frAeqindx = indxforcerows(indxforcerows >= lamindx.Aeqstart ...
                & indxforcerows <= lamindx.Aeqend) - numRowsAineq;
            
            rc = [ Aineq_orig(iAineqindex, fixforce);  Aeq_orig([iAeqindex csgrows], fixforce) ]' ...
                 * [ lambda.ineqlin(iAineqindex); lambda.eqlin([iAeqindex csgrows]) ] + c_orig(fixforce);
            A = [Aineq_orig(frAineqindx, fixforce);
                Aeq_orig(frAeqindx, fixforce)];
            
            yforce = zeros(nforcerows, 1);
            for i = 1: nforcerows
                tcol = A(i, :) ~= 0;
                if Rlbeqb(i) % in the case Rlb_i = b_i
                    ty = min( rc(tcol)./ A(i, tcol)' );
                else % in the case Rub_i = b_i  
                    ty = max( rc(tcol)./ A(i, tcol)' );
                end
                yforce(i) = ty;
                rc(tcol) = rc(tcol) - ty * A(i, tcol)';
            end
            
            lambda.ineqlin(frAineqindx) = - yforce(1: length(frAineqindx));
            lambda.eqlin(frAeqindx) = - yforce(length(frAineqindx)+1: end);
            
            lambda.lower(ffixedlower) = rc( 1:length(ffixedlower) );
            lambda.upper(ffixedupper) = -rc( length(ffixedlower)+1 : end ); 
        end %if restoreSol(4).exist
        
        % Case: row singeltons
        if restoreSol(5).exist && ~isempty(restoreSol(5).lorow)
            result = -c_orig - Aineq_orig'*lambda.ineqlin - Aeq_orig'*lambda.eqlin ...
                + lambda.lower - lambda.upper;
            sgrows = restoreSol(5).lorow;
            sgcols = restoreSol(5).upcol;
            lambda.eqlin(sgrows) = (result(sgcols))./diag(Aeq_orig(sgrows, sgcols));
        end
        
        % Case: fixed variables due to equal bounds
        if restoreSol(6).exist
            fixedlower = restoreSol(6).lorow;
            fixedupper = restoreSol(6).upcol;
            cplusAlambda = c_orig + Aineq_orig'*lambda.ineqlin + Aeq_orig'*lambda.eqlin;
            lambda.lower(fixedlower) =   max( cplusAlambda(fixedlower), 0);
            lambda.upper(fixedupper) =   max( -cplusAlambda(fixedupper), 0);
            % modify any negative assignment to lambda.lower and lambda.upper
            lambda.upper(fixedlower) = - min( cplusAlambda(fixedlower), 0);
            lambda.lower(fixedupper) =   max( cplusAlambda(fixedupper), 0);
        end
     
    end % if caseflag == 1
    
    if ( diagnostic_level >= 4 )
        disp( sprintf('  The norm of the dual feasibility is norm(-A''*y + w - z + c) = %8.2e.', ...
            norm(Aineq_orig'* lambda.ineqlin + Aeq_orig'*lambda.eqlin + lambda.upper -lambda.lower + c_orig)));
    end
    
end % if computeLambda

% Restore the basis index to the original system if required
% rbasicVarIdx: basic variable indices for the reduced problem after simplexpresolve.
% rnonbasicVarIdx: nonbasic variable indices for the reduced problem after simplexpresolve.
basicVarIdx = [];
nonbasicVarIdx = [];
delrows = [];

if restoreBasis
    if nnz(delrows1) >= 1
        MarkConstr(MarkConstr) = ~delrows1;
    end
    delrows = ~MarkConstr;
    varindx = find(MarkVars);
    basicVarIdx = varindx(rbasicVarIdx);
    nonbasicVarIdx = varindx(rnonbasicVarIdx);
    
    % Case: fix variables at bounds according to the objective coefficient
    % Add fixed variables to nonbasicVarIdx.
    if (restoreSol(8).exist)
        fixed = restoreSol(8).col; 
        nonbasicVarIdx = [nonbasicVarIdx; fixed];
    end
    
    % Case: generated zero column case
    % Add fixed variables to nonbasicVarIdx.
    if (restoreSol(7).exist)
        fixed = restoreSol(7).slvd_idx;
        nonbasicVarIdx = [nonbasicVarIdx; fixed];
    end

    % Case: free and implied free column singleton
    % Add column singleton variables to basicVarIdx and add the corresponding row.
    if (restoreSol(3).exist)
       solved = restoreSol(3).col; 
       basicVarIdx= [basicVarIdx; solved]; 
       delrows(restoreSol(3).row) = false;
    end
       
    % Case: row singletons
    % Add the solved variables by row singletons to basicVarIdx and add the
    % singleton row back.
    if (restoreSol(5).exist)
        solved = restoreSol(5).col;
        basicVarIdx = [basicVarIdx; solved(:)];
        delrows(restoreSol(5).row) = false;
    end
    
    % Case: zero column
    % Add the indices of fixed variables to nonbasicVarIdx.
    if (restoreSol(2).exist)
        fixed = restoreSol(2).slvd_idx;
        nonbasicVarIdx = [nonbasicVarIdx; fixed];
    end

    % Case: fixed variables by equal upper and lower bounds
    % Add the index of fixed variables to nonbasicVarIdx.
    if (restoreSol(6).exist)
        fixed = find(restoreSol(6).slvd_idx);
        nonbasicVarIdx = [nonbasicVarIdx; fixed];
    end         
end % if restorebasicVarIdx
