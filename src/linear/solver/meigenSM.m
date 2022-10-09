function [P, L, varargout] = meigenSM(K, M, bc, ne, Kold, R, psi0, s, options)
% [P, L] = MEIGENSM(K, M, bc, ne) finds the first ne eigenvalues to the 
% generalized eigenvalue problem KP = LMP with homogeneous boundary 
% conditions bc
%
% [P, L, R] = MEIGENSM(K, M, bc, ne) finds the first ne eigenvalues to the 
% generalized eigenvalue problem and returns R, the cholesky factorization 
% of the free part of K
%
% [P, L, deltas, B] = MEIGENSM(K, M, bc, ne, Kold, R, psi0, s, options)
% finds the first ne eigenvalues using CA. Kold is the stiffness matrix
% corresponding to the factorization R and the eigenmodes psi0. CA uses s
% basis vectors to compute the eigenvalues, and options is a struct
% containing information for the CA procedure. The basis vectors for the
% reduced order model are computed using CAeigs.
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
    % Solve the problem using cholesky factorization. MATLAB's eigs expects
    % the right hand side to be factorized, so we instead solve
    % M*V = (lambda)^(-1)K*V for the eigenparis whose eigenvalues have the
    % largest absolute value
    
    R = chol(Kff);
    [Pf, Linv] = eigs(Mff, R, ne, 'largestabs', ...
        'IsCholesky', true);
    L = diag(1./diag(Linv));
    
    % Normalize with respect to the mass matrix
    Pf = Pf./sqrt(dot(Pf, M*Pf));
    
    % Typically the user wants the cholesky factorization
    varargout = {R};
else
    % Solve the problem using a reduced order model if the proper
    % parameters are supplied
    psi0f = psi0(nf, :);
    dKff = Kff - Kold(nf, nf);
    
    [Pf, L, deltas, B] = CAeigs(Kff, Mff, ne, R, dKff, psi0f, s, options);
    varargout = {deltas, B};
end

% Sorting eigenvalues in ascending order
L = reshape(diag(L), ne, 1);
[Ld, I] = sort(L, 'ascend');
L = diag(Ld);
Pf = Pf(:, I);

% Making eigenvectors full
P = zeros(ndof, ne);
P(nf, :) = Pf;
end