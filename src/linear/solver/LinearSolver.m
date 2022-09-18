classdef LinearSolver < handle
    properties
        % Storing old factorizations and solutions
        Rold;                           % Old cholesky factorizations
        Kold;                           % Old stiffness matricies
        NSTEPS;                         % Number of steps for each equilibrium solve
        NSAVED;                         % Number of stiffness matricies to save
        
        % CA data
        ITSSF;                          % Number of iterations since factorization
        forceFactorization = 1;         % Boolean controlling if a factorization should be forced
        MAXITS;                         % Maximum number of opt steps before factorization
        NBASIS;                         % Number of basis vectors
        
        % Solver statistics such as number of factorizations and number of calls
        statistics;
    end
    
    methods
        function obj = LinearSolver(NBASIS, MAXITS, NSAVED, NSTEPS)
            obj.MAXITS = MAXITS;
            obj.ITSSF = MAXITS;
            
            obj.NBASIS  = NBASIS;
            obj.statistics = struct('CALLS', 0, 'FACTS', 0);
            
            obj.NSTEPS = NSTEPS;
            obj.NSAVED = NSAVED;
        end
        
        % Solves the equilibrium equations given by K, f, and bc
        function [x, f, fact] = solveq(obj, K, f, bc, n)
            % SOLVEQ solves the problem Ku = f with boundary conditions bc.
            % LinearSolver stores a number of stiffness matricies,
            % corresponding to the iteration number n.
            obj.statistics.CALLS = obj.statistics.CALLS + 1;
            
            % Default n is the last one
            if nargin < 5
                n = length(obj.Kold);
            elseif nargin < 4
                errorStruct.message = 'Too few arguments';
                error(errorStruct);
            end
            
            % Initializing the dofs that correspond to known displacements
            np = bc(:, 1);
            xp = bc(:, 2);
            
            ndof = length(f);
            nf = 1:ndof;
            nf(np) = [];
            
            % Initializes the stiffness matrix and extracts the submatrices
            % corresponding to solving the system
            % [A B][u] = [g]    u unknown,  g known
            % [C D][v] = [h]    v known,    h unknown
            Kff = K(nf, nf);
            Kfp = K(nf, np);
            Kpf = K(np, nf);
            Kpp = K(np, np);
            
            % Solving the reduced system Au = g - Bv
            ff = f(nf);
            b = ff - Kfp*xp;
            
            if obj.forceFactorization || ...
                    obj.ITSSF >= obj.MAXITS || ...
                    n <= obj.NSTEPS - obj.NSAVED || ...
                    n > length(obj.Kold)
                obj.ITSSF = 0;
                obj.statistics.FACTS = obj.statistics.FACTS + 1;
                fact = 1;
                
                % Try to do a cholesky, if the matrix is neg def do lu
                % instead.
                try
                    R = chol(Kff);
                    if n > obj.NSTEPS - obj.NSAVED
                        obj.Kold{n} = Kff;
                        obj.Rold{n} = R;
                    end
                    xf = R\(R'\b);
                catch
                    % Does not really work yet as i cannot save both U and
                    % L, so this just makes the method fail if the
                    % factorization is to be reused nest step.
                    
                    % maybe using a struct containing information of the
                    % factorization used could solve it.
                    [L, U] = lu(Kff);
                    xf = U\(L\b);
                end
                
            else
                % Reusing previous factorization of the submatrix A
                fact = 0;
                R = obj.Rold{n};
                dK = Kff - obj.Kold{n};
                
                % Generating basis vectors
                B = CASBON(R, dK, Kff, b, obj.NBASIS);
                V = B{end};
                
                % Projecting b onto the basis for the solution
                z = V'*b;
                xf = V*z;
            end
            % Compuing the reaction forces
            fp = Kpf*xf + Kpp*xp;
            f(np) = fp;
            
            % Assembling final solution vector
            x = zeros(ndof, 1);
            x(np) = xp;
            x(nf) = xf;
        end
        
        function flush(obj)
            obj.ITSSF = obj.ITSSF + 1;
            obj.statistics.CALLS = 0;
            obj.statistics.FACTS = 0;
        end
        
        function printStats(obj)
            fprintf('Linear Solver stats\n')
            fprintf(repmat('-', 1, 35));
            fprintf('\n');
            
            data = {obj.NBASIS, obj.MAXITS, obj.NSAVED, obj.NSTEPS};
            names = {'Basis vectors', 'Max its', 'Num Saved', 'Num Steps'};
            form = '%-18s %10i\n';
            for k = 1:numel(data)
                fprintf(form, names{k}, data{k});
            end
            
            fprintf(repmat('-', 1, 35));
            fprintf('\n');
        end
    end
end