classdef LinearSolver < handle
    properties
        % Storing old factorizations and solutions
        Rold;                           % Old cholesky factorizations
        Kold;                           % Old stiffness matricies
        nsteps;                         % Number of steps for each equilibrium solve
        nsaved = 2;                     % Number of stiffness matricies to save
        
        % CA data
        itssf;                          % Number of iterations since factorization
        forceFactorization = 1;         % Boolean controlling if a factorization should be forced
        maxits;                         % Maximum number of opt steps before factorization
        nbasis;                         % Number of basis vectors
        
        % Solver statistics such as number of factorizations and number of calls
        statistics;                         
    end
    
    methods
        function obj = LinearSolver(maxits, nbasis, nsteps)
            obj.maxits = maxits;
            obj.itssf = maxits;
            
            obj.nbasis  = nbasis;
            obj.statistics = struct('calls', 0, ...
                                    'facts', 0);
                                
            obj.nsteps = nsteps;
        end
        
        % Solves the equilibrium equations given by K, f, and bc
        function [x, f] = solveq(obj, K, f, bc, n)
            % SOLVEQ solves the problem Ku = f with boundary conditions bc.
            % LinearSolver stores a number of stiffness matricies,
            % corresponding to the iteration number n.
            obj.statistics.calls = obj.statistics.calls + 1;
            
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
            
            % Solve the reduced system Au = g - Bv
            ff = f(nf);
            b = ff - Kfp*xp;
            
            % If direct solvers is used maxits is 1, otherwise CA is used
            if obj.forceFactorization || ...
                    obj.itssf >= obj.maxits || ...
                    n <= obj.nsteps - obj.nsaved || ...
                    n > length(obj.Kold)
                obj.itssf = 0;
                obj.statistics.facts = obj.statistics.facts + 1;
                
                % Try to do a cholesky, if the matrix is neg def do lu
                % instead.
                try
                    R = chol(Kff);
                    if n > obj.nsteps - obj.nsaved
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
                % Reusing previous factorization R of the submatrix A
                R = obj.Rold{n};
                dK = Kff - obj.Kold{n};
                
                % Generating basis vectors
                B = CASBON(R, dK, Kff, b, obj.nbasis);
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
            obj.statistics.calls = 0;
            obj.statistics.facts = 0;
            obj.itssf = obj.itssf + 1;
        end
        
    end
end