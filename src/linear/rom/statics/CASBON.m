function B = CASBON(R, dK, K, f, s)
% CASBON computes orthonormal basis vectors for static problem
%
% INPUT:
%       R:      Cholesky factorization of stiffness matrix          (n x n)
%       dK:     Change in stiffness matrix                          (n x n)
%       K:      Stiffness matrix                                    (n x m)
%       f:      Right hand side                                     (n x 1)    
%       s:      Number of basis vectors to generate                 1
% 
% OUTPUT:
%       V:      Basis vectors                                       (n x s)
%

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
B = {U, T, Rb, V};