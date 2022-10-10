function [V, B] = CAEEON(R, dK, K, M, psi0, vectors, s)
% The number of vectors to orthogonalize to
m = size(vectors, 2);

% Compute first basis vector
ui = R\(R'\(M*psi0));
ti = ui/sqrt(ui'*M*ui);
vi = ti;
for j = 1:m
    psij = vectors(:, j);
    vi = vi - (ti'*M*psij)*psij;
end
U = ui;
T = ti;
V = vi;

% Compute remaining basis vectors
for i = 2:s
    ui = -R\(R'\(dK*ti));
    
    % Normalize
    ti = ui/sqrt(ui'*M*ui);
    
    % Orthogonalize
    vi = ti;
    for j = 1:m
        psij = vectors(:, j);
        vi = vi - (ti'*M*psij)*psij;
    end
    
    U(:, i) = ui;
    T(:, i) = ti;
    V(:, i) = vi;
end

% Return all basis vectors, which is required for a consistent sensitivity
% analysis
B = {U, T};
end