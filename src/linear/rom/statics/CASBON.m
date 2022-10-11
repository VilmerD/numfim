function [V, U, T, R] = CASBON(A, b, Rold, Aold, s)
% V = CASBON(A, b, Rold, Aold, s) computes s orthonormal basis vectors using CA for
% the static problem A*x = b. R is the cholesky factorization of the old matrix Aold.
% 
% [V, B] = CASBON(A, b, Rold, Aold, s) computes the basis vectors and returns the 
% intermediate basis vectors required for a consistent sensitivity analysis.
% Compute change in A
dA = A - Aold;

% Compute first basis vector 
ui = Rold\(Rold'\b);
ti = ui/sqrt(ui'*A*ui);
U = ui;
T = ti;
R = ti;
V = ti;

% Compute remaining basis vectors
for i = 2:s
    ui = -Rold\(Rold'\(dA*ti));
    ti = ui/sqrt(ui'*A*ui);
    
    % Orthogonalize
    ri = ti;
    for j = 1:(i-1)
        vj = V(:, j);
        ri = ri - (ti'*A*vj)*vj;
    end
    
    % Normalize
    vi = ri/sqrt(ri'*A*ri);
    
    % Insert
    U(:, i) = ui;
    T(:, i) = ti;
    R(:, i) = ri;
    V(:, i) = vi;
end
end