function B = CAERS(R, dK, K, M, psi0, s)
% CAERS computes basis vectors for the reduced eigenproblem using Rayleigh
% shift

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
end