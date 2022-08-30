function [u, Q, R, varargout] = msolveq(K, f, bc, Kold, R, s)
% MSOLVEQ solves the linear system of equations with boundary conditions,
% with the option to use CA if R, Kold and s are supplied
%
% Input:
%       K:      Stiffness matrix                                    (n x n)
%       f:      Forces                                              (n x 1)
%       bc:     Boundary Conditions                                 (m x 2)
%       R:      Cholesky Factorization                              (n x n)
%       Kold:   Old stiffness matrix                                (n x n)
%       s:      Number of basis vectors to generate                       1
%
% Output:
%       u:          Displacements                                   (n x 1)
%       Q:          Forces                                          (n x 1)
%       R:          Cholesky factorization                          (n x n)
%       varargout:  Contains basis vectors of CA                          1
%
% NOTE: - If K, f and bc only are supplied the system is solved exactly.
%       - If R, Kold and s are supplied in addition to K, f, and bc the
%       system is solved using combined approximations with basis
%       orthonormalization (CASBON)
%

% If bc is empty just solve directly
if isempty(bc)
    R = chol(K);
    u = R\(R'\f);
    Q = f;
    return
end

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
else
    % Using CA
    dK = Kff - Kold(nf, nf);
    B = CASBON(R, dK, Kff, ff, s);
    U = B{end};
    uf = U*(U'*ff);
    varargout = {B};
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