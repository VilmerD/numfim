function [u, Q, R, nits] = itersolveq(K, f, bc, R)
% ITERSOLVEQ solves the linear system of equations with boundary conditions,
% whith cholesky or pCG if symmetric preconditioner R is supplied
%
% Input:
%       K:      Stiffness matrix                                    (n x n)
%       f:      Forces                                              (n x 1)
%       bc:     Boundary Conditions                                 (m x 2)
%       R:      Cholesky Factorization                              (n x n)
%
% Output:
%       u:          Displacements                                   (n x 1)
%       Q:          Forces                                          (n x 1)
%       R:          Cholesky factorizati9on                          (n x n)
%

% Split up the system into free and prescribed nodes
ndof = size(K, 1);
nf = (1:ndof)';
np = bc(:, 1);
nf(np) = [];

Kff = K(nf, nf);
Kfp = K(nf, np);
Kpp = K(np, np);
up = bc(:, 2);
ff = f(nf) - Kfp*up;

% Solving free system
if nargin == 3
    % Exactly
    R = chol(Kff);
    uf = R\(R'\(ff));
    nits = 0;
else
    % Using pCG
    [uf, nits] = pCG(@(x) Kff*x, ff, R, 1e-12);
end

% Reassembling the solution
u = zeros(ndof, 1);
u(np) = up;
u(nf) = uf;

f = zeros(ndof, 1);
fp = Kfp'*uf + Kpp*up;
f(np) = fp;
f(nf) = ff;

Q = f;
end