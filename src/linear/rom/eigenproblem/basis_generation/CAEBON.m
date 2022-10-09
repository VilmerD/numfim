function V = CAEBON(R, dK, K, M, psi0, s)
% CAEBON computes basis vectors for the reduced eigenproblem using CA

% Compute first basis vector 
ui = R\(R'\(M*psi0));
ti = ui/sqrt(ui'*M*ui);
vi = ti;
V = vi;

% Compute remaining basis vectors
for i = 2:s
    ui = -R\(R'\(dK*ti));
    ti = ti/sqrt(ti'*M*ti);
    
    % Orthogonalize
    vi = ti;
    for j = 1:(i-1)
        vj = V(:, j);
        vi = ti - (ui'*K*vj)*vj;
    end
    
    % Normalize
    V(:, i) = vi;
end