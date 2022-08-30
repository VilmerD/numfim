classdef LinearSolver < handle
    properties
        % Storing old factorizations and solutions
        Rold;
        Kold;
        nsteps;
        
        % CA data
        iterationsSinceFactorization;
        forceFactorization = 1;
        maxits;
        nbasis;
        
        % Solver statistics such as number of factorizations and number of
        % calls
        statistics;
    end
    
    methods
        function obj = LinearSolver(maxits, nbasis)
            obj.maxits = maxits;
            obj.iterationsSinceFactorization = maxits;
            
            obj.nbasis  = nbasis;
            obj.statistics = struct('ncalls', 0, ...
                                    'factorizations', 0);
        end
        
        % Solves the equilibrium equations given by K, f, and bc
        function [f, x] = solveq(obj, K, f, bc, n)
            obj.statistics.ncalls = obj.statistics.ncalls + 1;
            
            % Initializing the dofs that correspond to known displacements
            % (bck) and known forces(ddk)
            np = bc(:, 1);
            xp = bc(:, 2);
            
            ndof = length(f);
            nf = 1:ndof;
            nf(np) = [];
            
            % Initializes the stiffness matrix and extracts the submatrices
            % corresponding to solving the system
            % [A B][u] = [g]    u unknown,  g known
            % [C D][v] = [h]    v known,    h unknown
            [Kff, Kfp, Kpf, Kpp] = extractSubmatrices(K, np, nf);
            
            % Solve the reduced system Au = b
            ff = f(nf);
            up = xp;
            b = ff - Kfp*up;
            
            % If direct solvers is used maxits is 1, otherwise CA is used
            if obj.forceFactorization || ...
                    obj.iterationsSinceFactorization >= obj.maxits
                obj.iterationsSinceFactorization = 0;
                obj.statistics.factorizations = ...
                    obj.statistics.factorizations + 1;
                
                % Try to do a cholesky, if the matrix is neg def do lu
                % instead.
                try
                    R = chol(Kff);
                    if n > obj.nsteps - 3
                        obj.Kold{n} = Kff;
                        obj.Rold{n} = R;
                    end
                    uf = R\(R'\b);
                catch
                    % Does not really work yet as i cannot save both U and
                    % L, so this just makes the method fail if the
                    % factorization is to be reused nest step.
                    
                    % maybe using a struct containing information of the
                    % factorization used could solve it.
                    [L, U] = lu(Kff);
                    uf = U\(L\b);
                end
                
                % Reusing previous factorization R of the submatrix A
            else
                R = obj.Rold{n};
                dK = Kff - obj.Kold{n};
                % Generating basis vectors
                
                % If Knew = Kold, then dK = 0 and CA fails since the space
                % is spanned by only 1 vector.
                V = CA(R, Kff, dK, b, obj.nbasis);
                
                % Projecting b onto the basis for the solution
                z = V'*b;
                uf = V*z;
            end
            % Compuing the reaction forces
            fp = Kpf*uf + Kpp*up;
            f(np) = fp;
            
            % Assembling final solution vector
            x = zeros(ndof, 1);
            x(np) = up;
            x(nf) = uf;
        end
        
        function stats = getStats(obj)
            % Fetch statistics
            stats = obj.statistics;
            
            % Reset statistics
            obj.statistics.ncalls = 0;
            obj.statistics.factorizations = 0;
            obj.forceFactorization = 0;
            obj.iterationsSinceFactorization = ...
                obj.iterationsSinceFactorization + 1;
        end
        
    end
end