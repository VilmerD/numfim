classdef NonlinearSolver < handle
    % NONLINEARSOLVER solves nonlinear finite element problems
    
    properties
        % Model object
        model;
        
        % Linear solver
        linear_solver;
        
        % Solver options
        TOT_LOAD;
        LOAD_STEPS;
        RESTARTS_MAX;
        Z_STEPS;
        RMAX = 1e-3;
        RMIN = 1e-8;
        
        % Solver data
        uold;
        
        % Data/statistics
        datanames;
        statistics;
    end
    
    methods (Access = public)
        function obj = NonlinearSolver(model, solver, TOT_LOAD, LOAD_STEPS, ...
                RESTARTS_MAX, Z_STEPS)
            % Put model
            obj.model = model;
            obj.uold = zeros(model.ndof, LOAD_STEPS);
            
            % Initialize values
            obj.TOT_LOAD = TOT_LOAD;
            obj.RESTARTS_MAX = RESTARTS_MAX;
            obj.LOAD_STEPS = LOAD_STEPS;
            obj.Z_STEPS = Z_STEPS;
            
            % Use default linear solver if none specified
            if nargin < 1
                solver = @msolveq;
            end
            obj.linear_solver = solver;
            
            % Set up stats
            obj.datanames = {'NR_FLAG', 'N_INNER', 'N_LOAD_STEPS', ...
                'RFREE', 'RTOT', 'FACTS', 'START_AT_ZERO', 'ATTEMPT_Z_STEP'};
            emptydata = cell(1, numel(obj.datanames));
            obj.statistics = cell2struct(emptydata, obj.datanames, 2);
        end
        
        function [u, P, ef, es] = solve(obj, xnew, xold, Kold, drdxdx)
            % Give pointer to function for ease
            solveq_ptr = @(K, f, bc, n) obj.linear_solver.solveq(K, f, bc, n);
            
            % Initialize flags and counters
            PROBLEM_SOLVED = 0;
            START_AT_ZERO = 0;
            ATTEMPT_Z_STEP = 0;
            
            RESTARTS = 0;
            FORCED = 0;
            
            % Initialize stiffness matrix etc
            dx = xnew - xold;
            sfun = @(ed) obj.model.deform(ed);
            
            % Clear iteration data
            for k=1:length(obj.statistics)
                obj.statistics(1) = [];
            end
            
            % Restarts NRDC at steps u0
            while ~PROBLEM_SOLVED
                
                % Get initial guess
                if START_AT_ZERO
                    % If we start at zero initial guess should be zeros
                    u0 = zeros(size(obj.uold, 1), 1);
                    obj.linear_solver.forceFactorization = 1;
                    
                elseif ATTEMPT_Z_STEP
                    % Else if z-step start at last solution
                    u0 = obj.uold(:, end);
                    % Can we attempt z-stepping without factorization
                    if FORCED
                        obj.linear_solver.forceFactorization = 0;
                    end
                else
                    % Else jump, starting at last solution
                    u0 = obj.uold(:, end - RESTARTS);
                    if RESTARTS > 0 && ~obj.linear_solver.forceFactorization
                        FORCED = 1;
                        obj.linear_solver.forceFactorization = 1;
                    end
                    
                end
                
                % Use design changes to guess disp field
                if RESTARTS == 0 && nargin == 5
                    bc0 = [obj.model.nf, zeros(length(obj.model.nf), 1)];
                    u0dz = solveq_ptr(Kold, -drdxdx, bc0, obj.LOAD_STEPS);
                    u0 = u0+u0dz;
                end
                
                % Try to solve the problem
                if ~ATTEMPT_Z_STEP || START_AT_ZERO
                    % By stepping in u
                    % Initialize stiffness and residual
                    Kx = @(ef, es) obj.model.K(ef, es, xnew);
                    rx = @(ef, es) obj.model.fint(ef, es, xnew);
                    
                    % Attempt a solution starting at u0
                    [u, P, ef, es, NR_FLAG, N_INNER, N_LOAD_STEPS, ...
                        RFree, RTot, Facts] = ...
                        NRDC(Kx, rx, sfun, obj.TOT_LOAD, u0, obj.LOAD_STEPS, ...
                        solveq_ptr, obj.RMIN);
                    
                    % Update quantities
                    data = {NR_FLAG, N_INNER, N_LOAD_STEPS, RFree, ...
                        RTot, Facts, START_AT_ZERO, ATTEMPT_Z_STEP};
                    obj.saveStats(data);
                    
                else
                    % By stepping in z, starting at last solution u
                    xi = xold;
                    for i = 0:(obj.Z_STEPS-1)
                        % Increment z and get corresponding stiffness
                        xi = xi + dx/obj.Z_STEPS;
                        Kxi = @(ef, es) obj.model.K(ef, es, xi);
                        rxi = @(ef, es) obj.model.fint(ef, es, xi);
                        
                        % Adjust rtol for iteration number
%                         RTOL = obj.RMAX*(obj.RMIN/obj.RMAX)^(i/obj.Z_STEPS);
                        RTOL = obj.RMIN;
                        
                        % Attempt a solution
                        [u, P, ef, es, NR_FLAG, N_INNER, N_LOAD_STEPS, ...
                            RFree, RTot, Facts] = ...
                            NRDC(Kxi, rxi, sfun, obj.TOT_LOAD, u0, ...
                            obj.LOAD_STEPS, solveq_ptr, RTOL);
                        
                        % Update quantities
                        data = {NR_FLAG, N_INNER, N_LOAD_STEPS, RFree, ...
                            RTot, Facts, START_AT_ZERO, ATTEMPT_Z_STEP};
                        obj.saveStats(data);
                        
                        % Check if problem was solved, if not exit
                        if NR_FLAG ~= 0; break; end
                        u0 = u(:, end);
                    end
                    
                    % Check if stepping in z solved the issues
                    if NR_FLAG == 0
                        PROBLEM_SOLVED = 1;
                    end
                end
                
                % Check results
                if NR_FLAG == 0
                    % Solution found! Get out!
                    PROBLEM_SOLVED = 1;
                    
                else
                    % Restart with other initial guess
                    RESTARTS = RESTARTS + 1;
                    
                    if RESTARTS <= obj.RESTARTS_MAX
                        continue
                    elseif ~ATTEMPT_Z_STEP && obj.Z_STEPS > 1
                        % Try iterating through z instead
                        ATTEMPT_Z_STEP = 1;
                    elseif ~START_AT_ZERO
                        % If restarting fails, start at zero
                        START_AT_ZERO = 1;
                    else
                        break;
                    end
                end
            end
            
            % Add the last nonzero vectors of u into uold
            FIRST_NONZERO = find(normAlong(u, 2, 1) ~= 0, 1, 'first');
            obj.uold(:, FIRST_NONZERO:end) = u(:, FIRST_NONZERO:end);
            
            % If solution is still not found there is a problem...
            if ~PROBLEM_SOLVED
                stats = obj.statistics
                errorStruct.message = 'Couldnt solve problem';
                error(errorStruct);
            end
        end
        
        function saveStats(obj, stats)
            l = size(obj.statistics, 2);
            obj.statistics(l+1) = cell2struct(stats, obj.datanames, 2);
        end
        
        function printStats(obj)
            fprintf('Nonlinear Solver stats\n')
            fprintf(repmat('-', 1, 35));
            fprintf('\n');
            
            tl = obj.TOT_LOAD;
            perc_load = 100*tl(find(tl(:, 2) ~= 0, 1, 'first'), 2);
            data = {perc_load, obj.LOAD_STEPS, obj.RESTARTS_MAX, ...
                obj.Z_STEPS};
            names = {'Percentage load', 'Load Steps', 'Restarts max', 'Z Steps'};
            form = '%-18s %10i\n';
            for k = 1:numel(data)
                fprintf(form, names{k}, data{k});
            end
            
            fprintf(repmat('-', 1, 35));
            fprintf('\n');
            
            obj.linear_solver.printStats();
        end
    end
end