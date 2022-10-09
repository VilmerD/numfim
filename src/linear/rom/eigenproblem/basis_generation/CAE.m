function B = CAE(R, dK, K, M, psi0, s)

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