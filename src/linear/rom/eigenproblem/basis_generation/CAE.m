function B = CAE(R, dK, K, M, psi0, s)
% CA computes basis vectors for the reduced eigenproblem using CA
%
% INPUT:
%       R:      Cholesky factorization of stiffness matrix          (n x n)
%       dK:     Change in stiffness matrix                          (n x n)
%       K:      Stiffness matrix
%       M:      Mass matrix                                         (n x n)    
%       psi0:   Previous eigenvector                                (n x 1)
%       s:      Number of basis vectors to generate                 1
% 
% OUTPUT:
%       B:      Cell of basis vectors                                     2
%

% Compute first basis vector 
ui = R\(R'\(M*psi0));
vi = ui/sqrt(ui'*M*ui);
U = ui;
V = vi;

% Compute remaining basis vectors
for i = 2:s
    ui = -R\(R'\(dK*vi));
    vi = ui/sqrt(ui'*M*ui);
    
    U(:, i) = ui;
    V(:, i) = vi;
end

B = {U, V};