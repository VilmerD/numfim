classdef EigenSolver < handle
    properties
        % Solver variables
        ne;
        nf;
        np;
        ndof;
        edof;
        enod;
        nelm;
        
        % CA variables
        factfreq = 5;
        itfact = 11;
        zfact;
        alphtol = 10e-2;
        
        % Orthogonalization vars
        bool_orth = 1;
        type_orth = 'current';
        bool_shift = 0;
        
        % Tolerance / stop criteria vars
        nbmax = 6;
        btol = 1e-3;
        stol = 1e-3;
        
        % Data
        Kold;
        Rold;
        Pold;
        vecorth;
        
        datapoints = {'Pr', 'relDeltaMax'};
        data;
        
        
        % Logging
        logger;
    end
    
    methods
        function obj = EigenSolver(model, ne)
            % Problem and mesh data
            obj.nf = model.nf;
            obj.np = model.np;
            obj.edof = model.edof;
            obj.enod = model.enod;
            obj.ndof = model.ndof;
            obj.nelm = model.nelm;
            obj.ne = ne;
            
            obj.zfact = zerose(model.nelm, 1);
            
            % Data to be saved
            emptydata = cell(1, numel(obj.datapoints));
            obj.data = cell2struct(emptydata, obj.datapoints, 2);
            
            % Logger
            
        end
        
        
        function [P, L, dLdz] = eigenSM(obj, K, M, K0, M0, zk)
            % Solves gen. eigenproblem K*V = D*M*V
            Kff = K(obj.nf, obj.nf);
            Mff = M(obj.nf, obj.nf);
            
            % Determine if CA should be used at all
            alph = dangle(zk, obj.zfact, false);
            if alph < obj.alphtol && obj.itfact < obj.factfreq
                % Build and solve reduced model
                obj.itfact = obj.itfact+1;
                [P, L, dLdz] = obj.solveReduced(Kff, Mff, K0, M0, zk);
            else
                % Solve full model
                obj.itfact = 1;
                [P, L, dLdz] = obj.solveFull(Kff, Mff, K0, M0, zk);
                obj.zfact = zk;
            end
            
            % Sorting eigenvalues
            [Ld, I] = sort(diag(L), 'ascend');         % OBS ASCENDING ORDER
            L = diag(Ld);
            P = P(:, I);
            dLdz = dLdz(I, :);
        end
        
        function [Pf, L, dLdz] = solveReduced(obj, K, M, K0, M0, zk)
            % Solving the reduced problem for each eigenvector
            DK = K-obj.Kold;
            Pf = zeros(obj.ndof, obj.ne);
            Pfk = zeros(obj.ndof, 1);
            L = zeros(obj.ne);
            
            dLdz = zeros(obj.ne, length(zk));
            
%             Iref = logical((zk > 0.1).*(zk < 0.95));
            Prs = cell(obj.ne, 1);
            MRDs = cell(obj.ne, 1);
            for k = 1:obj.ne
                % Build first basis vector
                ui = cholsolve(obj.Rold, M*obj.Pold(:, k));
                ti = ui/sqrt(ui'*M*ui);
                vi = gsorth(ti, obj.vecorth(:, 1:(k-1)), M);
                V = vi;
                
                % Solve reduced problem
                Kr = V'*K*V;
                Mr = V'*M*V;
                [Pr, Lk] = eigs(Kr, Mr, 1, 'smallestabs', 'IsCholesky', false);
                Pk = V*Pr/sqrt(Pr'*Mr*Pr);      % Normalize wrt Mass matrix
                
                % Compute sensitivity of current eigenpair
                Pfk(obj.nf, 1) = Pk;    % full eigenmode must be used, ie including bcs
                dLkdzi = Deigen(Pfk, Lk, K0, M0, obj.edof);
                
                % Expand space until tolerance is met
                space_accepted = false;
                MRD = zeros(1, obj.nbmax);
                while ~space_accepted
                    % Add one more basis vector
                    ui = -cholsolve(obj.Rold, DK*ti);
                    ti = ui/sqrt(ui'*M*ui);
                    vi = gsorth(ti, obj.vecorth(:, 1:(k-1)), M);
                    V = [V vi];
                    
                    % Solve reduced problem
                    Kr = V'*K*V;
                    Mr = V'*M*V;
                    [Pr, Lk] = eigs(Kr, Mr, 1, 'smallestabs');
                    Pk = (V*Pr)/sqrt(Pr'*Mr*Pr);
                    
                    % Compute sensitivities
                    Pfk(obj.nf, 1) = Pk;
                    dLkdzip1 = Deigen(Pfk, Lk, K0, M0, obj.edof);
                    
                    % Check condition
                    relchange = abs((dLkdzip1 - dLkdzi)./dLkdzi);
                    change_measure = geomean(relchange);
                    space_accepted = ...
                        change_measure < obj.stol || ...
                        abs(Pr(end)) < obj.btol || ...
                        size(V, 2) >= obj.nbmax;
                    
                    % Update variables
                    dLkdzi = dLkdzip1;
                    MRD(size(V, 2)) = change_measure;
                end
                % Insert solution
                L(k, k) = Lk;
                Pf(:, k) = Pfk;
                dLdz(k, :) = dLkdzi;
                
                % Update data
                MRDs{k} = MRD;
                Prs{k} = Pr;
                
                % Log
                nb = size(V, 2);
                prq = abs(Pr(end)/Pr(1));
                mrdend = MRD(nb);
                fprintf('%2i | %2i %12.6f %12.6f\n', k, nb, mrdend, prq);
                
                % Update vectors used in orthogonalization
                if strcmpi(obj.type_orth, 'current')
                    obj.vecorth(:, k) = Pk;
                end
            end
            
            obj.data(length(obj.data)+1) = cell2struct({Prs, MRDs}, obj.datapoints, 2);
        end
        
        function [P, L, dLdz] = solveFull(obj, K, M, K0, M0, zk)
            % Solve full problem
            obj.Rold = chol(K);
            obj.Kold = K;
            [Pr, Linv] = eigs(M, obj.Rold, obj.ne, 'largestabs', 'IsCholesky', true);
            L = diag(1./diag(Linv));
            
            % Normalize wrt mass matrix
            for k = 1:obj.ne
                Prk = Pr(:, k);
                Pr(:, k) = Prk/sqrt(Prk'*M*Prk);
            end
            
            % Make full
            P(obj.nf, :) = Pr;
            P(obj.np, :) = 0;
            
            % Compute sensitivities
            dLdz = Deigen(P, L, K0, M0, obj.edof);
            
            % Update vectors used in orthogonalization
            if strcmpi(obj.type_orth, 'old')
                obj.vecorth = P;
            end
            obj.Pold = P(obj.nf, :);
        end
        
    end
    
end