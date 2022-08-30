function U = CAEBON(R, dK, K, M, psi0, s)
% CABON computes basis vectors for the reduced eigenproblem using CA
% and orthonormalizes (Grahm-Schmidt) them
%
% INPUT:
%       R:      Cholesky factorization of stiffness matrix          (n x n)
%       dK:     Change in stiffness matrix                          (n x n)
%       K:      Stiffness matrix                                    (n x n)
%       M:      Mass matrix                                         (n x n)    
%       psi0:   Previous eigenvectors                               (n x m)
%       s:      Number of basis vectors to generate                       1
% 
% OUTPUT:
%       U:      Basis vectors                                       (n x s)
%

% Compute first basis vector 
uhat1 = R\(R'\(M*psi0));
u1 = matrixNormalize(uhat1, M);
U = u1;

% Compute remaining basis vectors
for i = 2:s
    uhati = -R\(R'\(dK*U(:, i-1)));
    
    % Orthogonalize
    for j = 1:(i-1)
        uj = U(:, j);
        uhati = uhati - (uhati'*K*uj)*uj;
    end
    
    % Normalize
    U(:, i) = matrixNormalize(uhati, K);
end