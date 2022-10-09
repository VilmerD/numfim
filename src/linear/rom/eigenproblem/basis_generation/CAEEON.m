function B = CAEEON(R, dK, K, M, psim, s)
% CAEEON computes basis vectors for the reduced eigenproblem using CA
% and orthogonalizes (Grahm-Schmidt) them to the eigenvectors in psi0

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