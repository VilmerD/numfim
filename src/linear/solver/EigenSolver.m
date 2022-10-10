classdef EigenSolver < handle
    properties
        % Solver variables
        ne;
        nf;
        np;
        
        % CA variables
        xfact;
        
        % Orthogonalization vars
        bool_orth = 1;
        type_orth = 'current';
        
        % Tolerance / stop criteria
        factfreq = 25;                  % Factorization frequency
        itfact = inf;                   % Iterations since factorization
        alphtol = 5e-2;                 % Tolerance in design changes
        stol = 1e-2;                	% Tolerance in sensitivity change
        
        nbmin = 2;                      % Minimum number of basis vectors
        nbmax = 8;                      % Maximum number of basis vectors
            
        % Data from factorization
        Kold;
        Rold;
        Pold;
        vecorth;
        
        % Previous solution
        xprev;
        Pprev;
        Lprev;
        dLdzprev;
        
        % Data saving
        data;
        
        % Logging
        logger;
    end
    
    methods
        function obj = EigenSolver(nf, np, ne)
            % Problem and mesh data
            obj.nf = nf;
            obj.np = np;
            obj.ne = ne;
                        
            % Data to be saved
            obj.data = DataSaver({'zchange', 'schange', 'bchange'});
        end
        
        
        function [P, L, dLdz] = eigenSM(obj, K, M, xk, sensfun)
            % Solves gen. eigenproblem K*V = D*M*V
            % Check if problem has been solved already
            Kff = K(obj.nf, obj.nf); Mff = M(obj.nf, obj.nf);
            
            % Compute the angle to the previous factorization
            if isempty(obj.xfact)
                s = 1; obj.xfact = xk;
            else
                s = sang(xk, obj.xfact);
            end
            
            % If the angle is below the tolerance and the factorization is
            % fresh enough use CA, else solve full problem
            if s < obj.alphtol && obj.itfact < obj.factfreq
                % Build and solve reduced model
                [P, L, dLdz, bchng, schng] = ...
                    obj.solveReduced(Kff, Mff, sensfun);
                obj.itfact = obj.itfact + 1;
            else
                % Solve full model
                [P, L, dLdz] = obj.solveFull(Kff, Mff, sensfun);
                bchng = []; schng = [];
                obj.xfact = xk; obj.itfact = 1;
            end
            
            % Saving data
            obj.xprev = xk;
            obj.data.saveData({s, bchng, schng});
        end
        
        function [Pf, L, dLdz, bchng, schng] = solveReduced(obj, K, M, ...
                sensfun)
            % Solving the reduced problem for each eigenvector
            DK = K-obj.Kold;
            Pf([obj.nf; obj.np], 1:obj.ne) = 0;
            Pfk = Pf(:, 1);
            L = zeros(obj.ne);
            dLdz = zeros(obj.ne, numel(obj.xfact));
            
            bchng = cell(obj.ne, 1); schng = cell(obj.ne, 1);
            for k = 1:obj.ne
                % Build first basis vector
                ui = cholsolve(obj.Rold, M*obj.Pold(:, k));
                ti = ui/sqrt(ui'*M*ui);
                vi = ti;
                if ~strcmpi(obj.type_orth, 'none')
                    vi = gsorth(vi, obj.vecorth(:, 1:(k-1)), M);
                end
                V = vi;
                
                % Solve reduced problem
                Kr = V'*K*V; Mr = V'*M*V;
                [Pr, Lk] = eigs(Kr, Mr, 1, 'smallestabs', ...
                    'IsCholesky', false);
                Pk = V*Pr/sqrt(Pr'*Mr*Pr);      % Full approximation
                
                % Compute sensitivity of current eigenpair
                Pfk(obj.nf, 1) = Pk;  
                dLkdzi = sensfun(Pfk, Lk);
                
                % Expand space until tolerance is met
                nb = 1;
                space_accepted = false;
                schange = zeros(1, obj.nbmax);
                while nb < obj.nbmin || (~space_accepted && nb < obj.nbmax)
                    % Add one more basis vector
                    nb = nb + 1;
                    ui = -cholsolve(obj.Rold, DK*ti);
                    ti = ui/sqrt(ui'*M*ui);
                    vi = ti;
                    if ~strcmpi(obj.type_orth, 'none')
                        vi = gsorth(vi, obj.vecorth(:, 1:(k-1)), M);
                    end
                    V = [V vi];
                    
                    % Solve reduced problem
                    Kr = V'*K*V; Mr = V'*M*V;
                    [Pr, Lk] = eigs(Kr, Mr, 1, 'smallestabs');
                    Pk = (V*Pr)/sqrt(Pr'*Mr*Pr);
                    
                    % Compute sensitivities
                    Pfk(obj.nf, 1) = Pk;
                    dLkdzip1 = sensfun(Pfk, Lk);
                    
                    % Compute relative change
                    abschange = abs(dLkdzip1 - dLkdzi);
                    relchange = abschange./abs(dLkdzip1);
                    
                    % Disregard changes in elements with low sensitivity
                    change = geomean(relchange, 2);
                    space_accepted = ...
                        change < obj.stol || ...
                        size(V, 2) >= obj.nbmax;
                    
                    % Update variables
                    dLkdzi = dLkdzip1;
                    schange(nb) = change;
                end
                % Insert solution
                L(k, k) = Lk;
                Pf(:, k) = Pfk;
                dLdz(k, :) = dLkdzi;
                
                % Update data
                schng{k} = schange;
                bchng{k} = Pr;
                
                % Log to user
                fprintf('%2i | %2i %12.6f %12.6f\n', k, nb, ...
                    schange(nb), abs(Pr(end)/Pr(1)));
                
                % Update vectors used in orthogonalization
                if strcmpi(obj.type_orth, 'current')
                    obj.vecorth(:, k) = Pk;
                end
            end
            [Pf, L, dLdz] = EigenSolver.sortEigpairs(Pf, L, dLdz);
        end
        
        function [Pf, L, dLdz] = solveFull(obj, K, M, sensfun)
            % Solve full problem
            % Compute and store cholesky factorization of K
            obj.Rold = chol(K);
            obj.Kold = K;
            [P, Linv] = eigs(M, obj.Rold, obj.ne, ...
                'largestabs', 'IsCholesky', true);
            L = diag(1./diag(Linv));
            
            % Normalize wrt mass matrix
            for k = 1:obj.ne
                Prk = P(:, k);
                P(:, k) = Prk/sqrt(Prk'*M*Prk);
            end
            
            % Make full
            Pf(obj.nf, :) = P;
            Pf(obj.np, :) = 0;
            
            % Compute sensitivities
            dLdz = sensfun(Pf, L);
            
            % Sort eigenpairs
            [Pf, L, dLdz] = EigenSolver.sortEigpairs(Pf, L, dLdz);
            
            % Update vectors used in orthogonalization
            if strcmpi(obj.type_orth, 'old')
                obj.vecorth = Pf(obj.nf, :);
            end
            obj.Pold = Pf(obj.nf, :);
        end
        
    end
    
    methods (Static)
        function [P, L, dLdz] = sortEigpairs(P, L, dLdz)
            % Sorting eigenpairs in ascending order 
            [Ld, I] = sort(diag(L), 'ascend');         
            P = P(:, I);
            L = diag(Ld);
            dLdz = dLdz(I, :);
        end 
    end
    
end