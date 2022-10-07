function [P, L, varargout] = meigenSM(K, M, bc, ne, Kold, R, psi0, s, options)
% [P, L] MEIGENSM(K, M, bc, ne) finds the first ne eigenvalues to the 
% generalized eigenvalue problem KP = LMP with boundary conditions bc
%
% [P, L, R] MEIGENSM(K, M, bc, ne) finds the first ne eigenvalues to the 
% generalized eigenvalue problem KP = LMP with boundary conditions bc where
% R is the cholesky factorization of the free part of K
%
%

% Extracting free part of vibrational problem
ndof = size(K, 1);
nf = (1:ndof)';
np = bc(:, 1);
nf(np) = [];

Kff = K(nf, nf);
Mff = M(nf, nf);

% Solving generalized eigenvalue problem
if nargin == 4
    % Fully
    R = chol(Kff);
    [Pf, Linv] = eigs(Mff, R, ne, 'largestabs', ...
        'IsCholesky', true);
    L = diag(1./diag(Linv));
    for i = 1:ne
        Vfi = Pf(:, i);
        Pf(:, i) = Vfi/sqrt(Vfi'*Mff*Vfi);
    end
    varargout = {R};
else
    % Using reduced model
    psi0f = psi0(nf, :);
    dKff = Kff - Kold(nf, nf);
    
    [Pf, L, deltas, B] = CAeigs(Kff, Mff, ne, R, dKff, psi0f, s, options);
    varargout = {deltas, B};
end

% Sorting eigenvalues
L = reshape(diag(L), ne, 1);
[Dd, I] = sort(L, 'ascend');     % OBS
L = diag(Dd);
Pf = Pf(:, I);

% Making eigenvectors full
P = zeros(ndof, ne);
P(nf, :) = Pf;
end