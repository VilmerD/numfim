function [P, L, varargout] = meigenSM(A, B, bc, ne, Rold, Bold, Pold, s, options)
% [P, L] = MEIGENSM(A, B, bc, ne) finds the first ne eigenvalues to the 
% generalized eigenvalue problem (A+xB)v = 0 with homogeneous boundary 
% conditions bc
%
% [P, L, R] = MEIGENSM(A, B, bc, ne) finds the first ne eigenvalues to the 
% generalized eigenvalue problem and returns R, the cholesky factorization 
% of the free part of K
%
% [P, L, d, B] = MEIGENSM(A, B, bc, ne, Rold, Bold, Pold, s, options)
% finds the first ne eigenvalues using CA. Kold is the stiffness matrix
% corresponding to the factorization R and the eigenmodes psi0. CA uses s
% basis vectors to compute the eigenvalues, and options is a struct
% containing information for the CA procedure. The basis vectors for the
% reduced order model are computed using CAeigs.
%

% Extracting free part of vibrational problem
ndof = size(A, 1);
nf = (1:ndof)';
np = bc(:, 1);
nf(np) = [];

Aff = A(nf, nf);
Bff = B(nf, nf);

% Solving generalized eigenvalue problem
if nargin <= 5
    % Solve the problem using cholesky factorization. 
    if nargin < 5
        Rcurr = chol(Aff);
    else
        Rcurr = Rold;
    end
    [Pf, L] = eigs(Aff, Rcurr, ne, 'largestabs', ...
        'IsCholesky', true);
    
    % Normalize with respect to the mass matrix
    Pf = Pf./sqrt(dot(Pf, Bff*Pf));
    
    % Typically the user wants the cholesky factorization
    varargout = {Rcurr};
else
    % Solve the problem using a reduced order model if the proper
    [Pf, L, d, B] = CAeigs(Aff, Bff, ne, Rold, Bold(nf, nf), Pold(nf, :), s, options);
    varargout = {d, B};
end

% Sorting eigenvalues in ascending order
L = reshape(diag(L), ne, 1);
[Ld, I] = sort(L, 'descend');
L = diag(Ld);
Pf = Pf(:, I);

% Making eigenvectors full
P = zeros(ndof, ne);
P(nf, :) = Pf;
end