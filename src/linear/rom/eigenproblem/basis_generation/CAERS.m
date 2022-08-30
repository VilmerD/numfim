function B = CAERS(R, dK, K, M, psi0, s)
% CAS computes basis vectors for the reduced eigenproblem using Rayleigh -
% shifted CA
%
% INPUT:
%       R:      Cholesky factorization of stiffness matrix          (n x n)
%       dK:     Change in stiffness matrix                          (n x n)
%       K:      Stiffness matrix                                    (n x n)
%       M:      Mass matrix                                         (n x n)    
%       psi0:   Previous eigenvector                                (n x 1)
%       s:      Number of basis vectors to generate                       1
% 
% OUTPUT:
%       B:      Basis vectors                                             2
%

% First basis vector
ui = R\(R'\(M*psi0));
vi = ui/sqrt(ui'*M*ui);
U = ui;
V = vi;

% Computing remaning basis vectors
for i = 2:s
    mu = (vi'*K*vi);                    % Rayleigh shift
    ui = -R\(R'\((dK*vi - mu*M*vi)));
    vi = ui/sqrt(ui'*M*ui);
    
    U(:, i) = ui;
    V(:, i) = vi;
end

B = {U, V};