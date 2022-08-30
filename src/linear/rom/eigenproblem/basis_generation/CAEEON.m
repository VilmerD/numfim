function B = CAEEON(R, dK, K, M, psim, s)
% CABON computes basis vectors for the reduced eigenproblem using CA
% and orthogonalizes (Grahm-Schmidt) them to the eigenvectors in psi0
%
% INPUT:
%       R:      Cholesky factorization of stiffness matrix          (n x n)
%       dK:     Change in stiffness matrix                          (n x n)
%       K:      Stiffness matrix                                    (n x n)
%       M:      Mass matrix                                         (n x n)    
%       psim:   Previous and lower order eigenvectors               (n x m)     
%       s:      Number of basis vectors to generate                       1 
% 
% OUTPUT:
%       B:      Basis vectors                                             3
%

psi0 = psim(:, 1);
psi = psim(:, 2:end);
mm1 = size(psi, 2);

% Compute basis vectors
% Compute first basis vector 
ui = R\(R'\(M*psi0));
ti = ui/sqrt(ui'*M*ui);
vi = ti;
for j = 1:mm1
    psij = psi(:, j);
    vi = vi - (ti'*M*psij)*psij;
end
U = ui;
T = ti;
V = vi;

% Compute remaining basis vectors
for i = 2:s
    ui = -R\(R'\(dK*ti));
    ti = ui/sqrt(ui'*M*ui);
    
    vi = ti;
    for j = 1:mm1
        psij = psi(:, j);
        vi = vi - (ti'*M*psij)*psij;
    end
    
    U(:, i) = ui;
    T(:, i) = ti;
    V(:, i) = vi;
end

B = {U, T, V};