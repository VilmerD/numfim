classdef EigenSolver < handle
    properties
        % Solver variables
        ne;
        nf;
        np;
        edof;
        nelm;
        
        % CA variables
        factfreq = 10;
        itfact = 11;
        zfact;
        alphtol = 5e-2;
        
        % Orthogonalization vars
        bool_orth = 1;
        type_orth = 'current';
        bool_shift = 0;
        
        % Tolerance / stop criteria vars
        nbmax = 4;
        btol = 1e-6;
        stol = 1e-4;
        
        % Data
        Kold;
        Rold;
        Pold;
        vecorth;
        
        datapoints = {'Pr', 'relDeltaMax'};
        data;
    end
    
    methods
        function obj = EigenSolver(model, ne, nb)
            obj.nf = model.nf;
            obj.np = model.np;
            obj.edof = model.edof;
            obj.nelm = model.nelm;
            obj.ne = ne;
            obj.nbmax = nb;
            
            obj.zfact = zerose(model.nelm, 1);
            
            emptydata = cell(1, numel(obj.datapoints));
            obj.data = cell2struct(emptydata, obj.datapoints, 2);
        end
        
        
        function [P, L, dLdz] = eigenSM(obj, K, M, K0, M0, zk)
            % Solves gen. eigenproblem K*V = D*M*V
            Kff = K(obj.nf, obj.nf);
            Mff = M(obj.nf, obj.nf);
            
            % Determine if CA should be used at all
            alph = dangle(zk, obj.zfact, false);
            if alph < obj.alphtol && obj.itfact <= obj.factfreq
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
            L = reshape(diag(L), obj.ne, 1);
            [L, I] = sort(L, 'ascend');         % OBS ASCENDING ORDER
            P = P(:, I);
            dLdz = dLdz(I, :);
        end
        
        function [P, L, dLdz] = solveReduced(obj, K, M, K0, M0, zk)
            % Solving the reduced problem for each eigenvector
            DK = K-obj.Kold;
            P([obj.nf; obj.np], 1) = 0;
            L = zeros(obj.ne);
            
            dLdz = zeros(obj.ne, length(zk));
            
            Iref = logical((zk > 0.1).*(zk < 0.95));
            Prs = cell(obj.ne, 1);
            MRDs = cell(obj.ne, 1);
            for k = 1:obj.ne
                % Build first basis vector
                u1 = cholsolve(obj.Rold, M*obj.Pold(:, k));
                t1 = u1/sqrt(u1'*M*u1);
                v1 = t1;
                for j = 1:k-1
                    voj = obj.vecorth(:, j);
                    v1 = v1 - (t1'*M*voj)*voj;
                end
                V = v1;
                tim1 = t1;
                
                % Solve reduced problem
                Kr = V'*K*V;
                Mr = V'*M*V;
                [Pr, Lr] = eigs(Kr, Mr, 1, 'smallestabs', 'IsCholesky', false);
                Pk = V*Pr/sqrt(Pr'*Mr*Pr);
                Lk = Lr;
                
                % Compute sensitivities
                P(obj.nf, k) = Pk;
                P(obj.np, k) = 0;
                dLdzi = Deigen(P, Lk, obj.edof, K0, M0);
                
                % Expand space until tolerance is met
                space_accepted = false;
                MRD = zeros(1, obj.nbmax);
                fprintf('%2i\n', k);
                while ~space_accepted
                    % Add one more basis vector
                    ui = -cholsolve(obj.Rold, DK*tim1);
                    ti = ui/sqrt(ui'*M*ui);
                    
                    vi = ti;
                    for j = 1:k-1
                        voj = obj.vecorth(:, j);
                        vi = vi - (ti'*M*voj)*voj;
                    end
                    V = [V vi];
                    tim1 = ti;
                    
                    % Solve reduced problem
                    Kr = V'*K*V;
                    Mr = V'*M*V;
                    [Pr, Lr] = eigs(Kr, Mr, 1, 'smallestabs');
                    Pk = V*Pr/sqrt(Pr'*Mr*Pr);
                    
                    % Compute sensitivities
                    P(obj.nf, k) = Pk;
                    P(obj.np, k) = 0;
                    L(k, k) = Lr;
                    dLdzip1 = Deigen(P, Lr, obj.edof, K0, M0);
                    
                    % Check condition
                    relchange = abs((dLdzip1(Iref)) - dLdzi(Iref))./dLdzi(Iref);
                    maxrelchange = max(relchange);
                    space_accepted = ...
                        maxrelchange <= obj.stol || ...
                        abs(Pr(end)) < obj.btol || ...
                        size(V, 2) >= obj.nbmax;
                    
                    % Update variables
                    dLdzi = dLdzip1;
                    MRD(size(V, 2)) = maxrelchange;
                    fprintf('%12.8f %12.8f\n', maxrelchange, abs(Pr(end)));
                end
                dLdz(k, :) = dLdzi;
                MRDs{k} = MRD;
                Prs{k} = Pr;
            
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
            obj.Pold = Pr;
            
            % Normalize wrt mass matrix
            for k = 1:obj.ne
                Prk = Pr(:, k);
                Pr(:, k) = Prk/sqrt(Prk'*M*Prk);
            end
            
            % Make full
            P(obj.nf, :) = Pr;
            P(obj.np, :) = 0;
            
            % Compute sensitivities
            dLdz = Deigen(P, L, obj.edof, K0, M0);
            
            % Update vectors used in orthogonalization
            if strcmpi(obj.type_orth, 'old')
                obj.vecorth = P;
            end
        end
        
    end
    
end

function firstCA()

end

function nCA()

end