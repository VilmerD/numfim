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
        Z_GRANULARITY;
        
        % Solver data
        uold;
        
        % Data/statistics
        statistics;
    end
    
    methods
        function obj = NonlinearSolver(model, solver, TOT_LOAD, LOAD_STEPS, ...
                RESTARTS_MAX, Z_GRANULARITY)
            % Put model
            obj.model = model;
            obj.uold = zeros(model.ndof, LOAD_STEPS);
            
            % Initialize values
            obj.TOT_LOAD = TOT_LOAD;
            obj.RESTARTS_MAX = RESTARTS_MAX;
            obj.LOAD_STEPS = LOAD_STEPS;
            obj.Z_GRANULARITY = Z_GRANULARITY;
            
            % Use default linear solver if none specified
            if nargin < 1
                solver = @msolveq;
            end
            obj.linear_solver = solver;
            
            % Set up stats
            obj.statistics;
        end
        
        function [u, P, ef, es] = solve(obj, znew, zold)
            % Give pointer to function for ease
            solveq_ptr = @(K, f, bc, n) obj.linear_solver.solveq(K, f, bc, n);
            
            % Initialize flags and counters
            PROBLEM_SOLVED = 0;
            START_AT_ZERO = 0;
            RESTARTS = 0;
            TOT_LOAD_STEPS = 0;
            N_INNER_TOT = 0;
            ATTEMPT_Z_STEP = 0;
            
            % Initialize stiffness matrix etc
            dz = znew - zold;
            Kz = @(ef, es) obj.model.K(ef, es, znew);
            rz = @(ef, es) obj.model.fint(ef, es, znew);
            sfun = @(ed) obj.model.deform(ed);
            
            % Restarts NRDC at steps u0
            while ~PROBLEM_SOLVED
                
                % Get initial guess
                if RESTARTS < obj.RESTARTS_MAX
                    u0n = obj.uold(:, end - RESTARTS);
                    
                elseif ATTEMPT_Z_STEP
                    u0n = obj.uold(:, end);
                    
                elseif START_AT_ZERO
                    u0n = zeros(size(obj.uold, 1), 1);
                    obj.linear_solver.forceFactorization = 1;
                end
                
                % Try to solve the problem
                if ~ATTEMPT_Z_STEP
                    % By stepping in u
                    [u, P, ef, es, NR_FLAG, N_INNER, N_LOAD_STEPS] = ...
                        NRDC(Kz, rz, sfun, obj.TOT_LOAD, u0n, obj.LOAD_STEPS, ...
                        solveq_ptr);
                    N_INNER_TOT = N_INNER_TOT + N_INNER;
                    TOT_LOAD_STEPS = TOT_LOAD_STEPS + N_LOAD_STEPS;
                else
                    % By stepping in z, starting at last solution u
                    dzi = dz*obj.Z_GRANULARITY;
                    ui = u0n;
                    for i = 1:1/obj.Z_GRANULARITY
                        % Increment z and get corresponding stiffness
                        zi = zold + i*dzi;
                        Kzi = @(ef, es) obj.model.K(ef, es, zi);
                        rzi = @(ef, es) obj.model.fint(ef, es, zi);
                        
                        % Attempt a solution
                        [u, P, ef, es, NR_FLAG, N_INNER, N_LOAD_STEPS] = ...
                            NRDC(Kzi, rzi, sfun, obj.TOT_LOAD, ui, ...
                            obj.LOAD_STEPS, solveq_ptr);
                        N_INNER_TOT = N_INNER_TOT + N_INNER;
                        TOT_LOAD_STEPS = TOT_LOAD_STEPS + N_LOAD_STEPS;
                        
                        % Check if problem was solved, if not exit
                        if NR_FLAG ~= 0; break; end
                        ui = u(:, end);
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
                    
                elseif NR_FLAG == 1 && RESTARTS < obj.RESTARTS_MAX
                    % Restart with other initial guess
                    RESTARTS = RESTARTS + 1;
                    
                elseif NR_FLAG == 1 && ~ATTEMPT_Z_STEP && obj.Z_GRANULARITY < 1
                    % Try iterating through z instead
                    ATTEMPT_Z_STEP = 1;
                    
                elseif NR_FLAG == 1 && RESTARTS == obj.RESTARTS_MAX
                    % If restarting fails, start at zero
                    START_AT_ZERO = 1;
                    
                else
                    % Started from zero, and coulnd't solve problem. Need
                    % to exit the solver.
                    break;
                end
                
            end
            
            % Save some stats
            obj.statistics = struct(...
                'N_INNER_TOT', N_INNER_TOT, ...
                'N_RESTARTS', RESTARTS, ...
                'TOT_LOAD_STEPS', TOT_LOAD_STEPS, ...
                'ATTEMPT_Z_STEP', ATTEMPT_Z_STEP, ...
                'START_AT_ZERO', START_AT_ZERO);
            
            % Add the last nonzero vectors of u into uold
            FIRST_NONZERO = find(normAlong(u, 2, 1) ~= 0, 1, 'first');
            obj.uold(:, FIRST_NONZERO:end) = u(:, FIRST_NONZERO:end);
            
            % If solution is still not found there is a problem...
            if ~PROBLEM_SOLVED
                errorStruct.message = 'Couldnt solve problem';
                error(errorStruct);
            end
        end
    end
end