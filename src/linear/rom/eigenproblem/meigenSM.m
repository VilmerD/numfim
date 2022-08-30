function [V, D, varargout] = meigenSM(K, M, bc, ne, Kold, R, psi0, s, options)
% MEIGENSM solves the generalized eigenvalue problem with boundary
% conditions, and uses CA if right data is supplied
%
% Input:
%       K:      Stiffness matrix                                    (n x n)
%       M:      Mass matrix                                         (n x n)
%       bc:     Boundary Conditions                                 (m x 2)
%       k:      Number of eigenvalues to compute                          1
%       R:      Cholesky Factorization                              (n x n)
%       Kold:   Old stiffness matrix                                (n x n)
%       psi0:   Old eigenmodes                                      (n x k)
%       s:      Number of basis vectors to generate                       1
%       feig:   Basis generation function
%
% Output:
%       V:          Eigenmodes                                      (n x k)
%       D:          Eigenvalues                                     (k x 1)
%       varargout:  Contains basis vectors of CA                          1
%
% NOTE: - If K, M, bc and k only are supplied the system is solved exactly.
%       - If R, Kold, psi0, s and feig are supplied in addition to K, M, bc
%       and k the system is solved using combined approximations method
%       feig
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
    [Vf, Drec] = eigs(Mff, R, ne, 'largestabs', ...
        'IsCholesky', true);
    D = diag(1./diag(Drec));
    for i = 1:ne
        Vfi = Vf(:, i);
        Vf(:, i) = Vfi/sqrt(Vfi'*Mff*Vfi);
    end
    varargout = {R};
else
    % Using reduced model
    psi0f = psi0(nf, :);
    dKff = Kff - Kold(nf, nf);
    
    [Vf, D, deltas, B] = CAeigs(Kff, Mff, ne, R, dKff, psi0f, s, options);
    varargout = {deltas, B};
end
% Sorting eigenvalues
D = reshape(diag(D), ne, 1);
[D, I] = sort(D, 'ascend');     % OBS
Vf = Vf(:, I);

% Making eigenvectors full
V = zeros(ndof, ne);
V(nf, :) = Vf;
end