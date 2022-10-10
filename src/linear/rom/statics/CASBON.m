function [V, B] = CASBON(R, dK, K, f, s)
% V = CASBON(R, dK, K, f, s) computes s orthonormal basis vectors using CA for
% the static problem K*u = f. R is the cholesky factorization of the old
% stiffness matrix Kold = K + dK.
% 
% [V, B] = CASBON(R, dK, K, f, s) computes the basis vectors and returns the 
% intermediate basis vectors required for a consistent sensitivity analysis.

% Compute first basis vector 
ui = R\(R'\f);
ti = ui/sqrt(ui'*K*ui);
U = ui;
T = ti;
Rb = ti;
V = ti;

% Compute remaining basis vectors
for i = 2:s
    ui = -R\(R'\(dK*ti));
    ti = ui/sqrt(ui'*K*ui);
    
    % Orthogonalize
    ri = ti;
    for j = 1:(i-1)
        vj = V(:, j);
        ri = ri - (ti'*K*vj)*vj;
    end
    
    % Normalize
    vi = ri/sqrt(ri'*K*ri);
    
    % Insert
    U(:, i) = ui;
    T(:, i) = ti;
    Rb(:, i) = ri;
    V(:, i) = vi;
end

B = {U, T, Rb};
end