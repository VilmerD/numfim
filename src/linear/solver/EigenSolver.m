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
        type_orth = 'old';
        bool_shift = 0;
        
        % Tolerance / stop criteria
        factfreq = 5;           % Factorization frequency
        itfact = inf;           % Iterations since factorization
        nbmax = 4;              % Maximum number of basis vectors
        alphtol = 5e-2;         % Tolerance in design changes
        btol = 1e-2;            % Tolerance in y-components
        stol = 1e-2;            % Tolerance in sensitivity change
        
        % Data
        Kold;
        Rold;
        Pold;
        vecorth;
        
        xprev;
        Pprev;
        Lprev;
        dLdzprev;
        
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
            if ~isempty(obj.xprev) && norm(xk - obj.xprev) == 0
                P = obj.Pprev;
                L = obj.Lprev;
                dLdz = obj.dLdzprev;
                return;
            end
            Kff = K(obj.nf, obj.nf); Mff = M(obj.nf, obj.nf);
            
            % Determine if CA should be used at all
            if isempty(obj.xfact); s = 1; obj.xfact = xk;
            else; s = sang(xk, obj.xfact);
            end
            
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
            obj.Pprev = P;
            obj.Lprev = L;
            obj.dLdzprev = dLdz;
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
                
                if ~strcmpi(obj.type_orth, 'none')
                    vi = gsorth(ti, obj.vecorth(:, 1:(k-1)), M);
                else; vi = ti;
                end
                V = vi;
                
                % Solve reduced problem
                Kr = V'*K*V; Mr = V'*M*V;
                [Pr, Lk] = eigs(Kr, Mr, 1, 'smallestabs', 'IsCholesky', false);
                Pk = V*Pr/sqrt(Pr'*Mr*Pr);      % Normalize wrt Mass matrix
                
                % Compute sensitivity of current eigenpair
                Pfk(obj.nf, 1) = Pk;    % full eigenmode must be used, ie including bcs
                dLkdzi = sensfun(Pfk, Lk);
                
                % Expand space until tolerance is met
                nb = 1;
                space_accepted = false;
                schange = zeros(1, obj.nbmax);
                while ~space_accepted
                    % Add one more basis vector
                    nb = nb + 1;
                    ui = -cholsolve(obj.Rold, DK*ti);
                    ti = ui/sqrt(ui'*M*ui);
                    if ~strcmpi(obj.type_orth, 'none')
                        vi = gsorth(ti, obj.vecorth(:, 1:(k-1)), M);
                    else; vi = ti;
                    end
                    V = [V vi];
                    
                    % Solve reduced problem
                    Kr = V'*K*V; Mr = V'*M*V;
                    [Pr, Lk] = eigs(Kr, Mr, 1, 'smallestabs');
                    Pk = (V*Pr)/sqrt(Pr'*Mr*Pr);
                    
                    % Compute sensitivities
                    Pfk(obj.nf, 1) = Pk;
                    dLkdzip1 = sensfun(Pfk, Lk);
                    
                    % Check condition
                    relchange = abs((dLkdzip1 - dLkdzi)./dLkdzi);
                    change = geomean(relchange, 2);
                    space_accepted = ...
                        change < obj.stol || ...
                        abs(Pr(end)) < obj.btol || ...
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
                schng{k} = schange; bchng{k} = Pr;
                
                % Log
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