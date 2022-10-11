function [V, U, T] = CAEEON(A, B, Rold, Bold, Pold, s, VO)
% The number of vectors to orthogonalize to
m = size(VO, 2);
dB = B - Bold;

% Compute first basis vector
ui = Rold\(Rold'\(A*Pold));
ti = ui/sqrt(abs(ui'*B*ui));
vi = ti;
for j = 1:m
    vj = VO(:, j);
    vi = vi - (ti'*A*vj)*vj;
end
U = ui;
T = ti;
V = vi;

% Compute remaining basis vectors
for i = 2:s
    ui = -Rold\(Rold'\(dB*ti));
    
    % Normalize
    ti = ui/sqrt(abs(ui'*B*ui));
    
    % Orthogonalize
    vi = ti;
    for j = 1:m
        vj = VO(:, j);
        vi = vi - (ti'*A*vj)*vj;
    end
    
    U(:, i) = ui;
    T(:, i) = ti;
    V(:, i) = vi;
end
end