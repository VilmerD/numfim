function [V, D, deltas, B] = CAeigs(K, M, n, R, dK, psi0, s, options)
% CAeigs solves the generalized eigenproblem using a reduced order model
%
% INPUT:
%       K:      Stiffness matrix                                    (n x n)
%       M:      Mass matrix                                         (n x n)

% OUTPUT:
%       V:      Eigenvectors                                        (n x k)
%       D:      Eigenvalues                                         (k x k)

nelm = size(K, 1);
V = zeros(nelm, n);
D = zeros(n, n);
B = cell(n, 2);
deltas = -1*ones(n, 1);
for k = 1:n
    % Choose basis generation method
    vectors = [];
    switch options.orthotype
        case 'current'
            vectors = V(:, 1:(k-1));
        otherwise
            if ~isempty(options.orthovecs)
                vectors = options.orthovecs(:, 1:(k - 1));
            end
    end
    Bk = generateBasis(R, dK, K, M, psi0(:, k), vectors, ...
        options.shift, s);
    
    Vk = Bk{end};
    
    % Compute reduced model
    KR = Vk'*K*Vk;
    MR = Vk'*M*Vk;
    
    [yk, dk] = eigs(KR, MR, 1, 'smallestabs');  % Solve reduced model
    
    yk = yk/sqrt(yk'*MR*yk);                    % Normalize wrt mass matrix
    psik = Vk*yk;                               % Full solution
    
    V(:, k) = psik;                             % Insert solution
    D(k, k) = dk;
    B{k, 1} = Bk;
    B{k, 2} = yk;
    deltas(k) = norm(K*psik - dk*M*psik)/norm(K*psik);
end

end
function B = generateBasis(R, dK, K, M, psi0, vectors, shift, s)
orthogonalize = ~isempty(vectors);
mm1 = size(vectors, 2);

% Compute first basis vector
ui = R\(R'\(M*psi0));
ti = ui/sqrt(ui'*M*ui);
vi = ti;
for j = 1:mm1
    psij = vectors(:, j);
    vi = vi - (ti'*M*psij)*psij;
end
U = ui;
T = ti;
V = vi;

% Compute remaining basis vectors
for i = 2:s
    if shift
        mu = (vi'*K*vi);
        ui = -R\(R'\(dK*ti - mu*M*ti));
    else
        ui = -R\(R'\(dK*ti));
    end
    
    ti = ui/sqrt(ui'*M*ui);
    
    vi = ti;
    for j = 1:mm1
        psij = vectors(:, j);
        vi = vi - (ti'*M*psij)*psij;
    end
    
    U(:, i) = ui;
    T(:, i) = ti;
    V(:, i) = vi;
end

% Output args
if orthogonalize
    B = {U, T, V};
else
    B = {U, V};
end
end