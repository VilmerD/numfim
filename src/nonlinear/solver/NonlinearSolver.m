classdef NonlinearSolver < handle
    % NONLINEARSOLVER solves nonlinear finite element problems
    
    properties
        % Model object
        model;
        
        % Linear solver
        linear_solver;
        
        % Solver options
        RESTARTS_MAX;
        LOAD_STEPS;
        Z_GRANULARITY;
        
        % Solver data
        uold;
        
        % Data/statistics
        statistics;
    end
    
    methods
        function obj = NonlinearSolver(model, LOAD_STEPS, solver)
            % Put model
            obj.model = model;
            obj.uold = zeros(model.ndof, LOAD_STEPS);
            
            % Initialize values
            obj.RESTARTS_MAX = 3;
            obj.LOAD_STEPS = LOAD_STEPS;
            obj.Z_GRANULARITY = 0.25;
            
            % Use default linear solver if none specified
            if nargin < 1
                solver = @msolveq;
            end
            obj.linear_solver = solver;
            
            % Set up stats
            obj.statistics;
        end
        
        function [u, P, ef, es] = solve(obj, znew, zold, tot_load)
            % Give pointer to function for ease
            solveq_ptr = @(K, f, bc, n) obj.linear_solver.solveq(K, f, bc, n);
            
            % Initialize flags and counters
            PROBLEM_SOLVED = 0;
            START_AT_ZERO = 0;
            RESTARTS = 0;
            N_INNER_TOT = 0;
            
            % Initialize stiffness matrix etc
            dz = znew - zold;
            Kz = @(ef, es) obj.model.K(ef, es, znew);
            rz = @(ef, es) obj.model.fint(ef, es, znew);
            sfun = @(ed) obj.model.deform(ed);
            
            % Restarts NRDC at steps u0
            while ~PROBLEM_SOLVED
                
                % Get initial guess
                if ~START_AT_ZERO && RESTARTS <= obj.RESTARTS_MAX
                    u0n = obj.uold(:, end - RESTARTS);
                else
                    START_AT_ZERO = 1;
                    u0n = zeros(size(obj.uold, 1), 1);
                    obj.linear_solver.forceFactorization = 1;
                end
                
                % Attempt a solution at u0n
                [u, P, ef, es, NR_FLAG, N_INNER] = ...
                    NRDC(Kz, rz, sfun, tot_load, u0n, obj.LOAD_STEPS, ...
                    solveq_ptr);
                N_INNER_TOT = N_INNER_TOT + N_INNER;
                
                % Check results
                if NR_FLAG == 0
                    % Solution found! Get out!
                    PROBLEM_SOLVED = 1;
                    
                elseif NR_FLAG == 1 && RESTARTS < obj.RESTARTS_MAX
                    % Restart with other initial guess
                    RESTARTS = RESTARTS + 1;
                    
                elseif NR_FLAG == 1 && RESTARTS == obj.RESTARTS_MAX
                    % Try iterating through z instead
                    dzi = dz*obj.Z_GRANULARITY;
                    ui = obj.uold(:, end);
                    for i = 1:1/obj.Z_GRANULARITY
                        % Increment z and get corresponding stiffness
                        zi = zold + i*dzi;
                        Kzi = @(ef, es) obj.model.K(ef, es, zi);
                        rzi = @(ef, es) obj.model.fint(ef, es, zi);
                        
                        % Attempt a solution
                        [u, P, ef, es, NR_FLAG, N_INNER] = NRDC(Kzi, rzi, sfun, ...
                            tot_load, ui, obj.LOAD_STEPS, solveq_ptr);
                        N_INNER_TOT = N_INNER_TOT + N_INNER;
                        
                        % Check if problem was solved, if not exit
                        if NR_FLAG ~= 0; break; end
                        ui = u(:, end);
                    end
                    
                    % Check if stepping in z solved the issues
                    if NR_FLAG == 0
                        PROBLEM_SOLVED = 1;
                    else
                        RESTARTS = RESTARTS + 1;
                    end
                    
                elseif START_AT_ZERO
                    % Started from zero, and coulnd't solve problem. Need
                    % to exit the solver.
                    break;
                end
                
            end
            
            obj.statistics = struct(...
                'N_INNER_TOT', N_INNER_TOT, ...
                'N_RESTARTS', RESTARTS, ...
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