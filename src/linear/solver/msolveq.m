function [u, Q, R, varargout] = msolveq(K, f, bc, Kold, R, s)
% [u, Q, R] = MSOLVEQ(K, f, bc) solves the linear system of equations 
% K*u = f with boundary conditions bc. Q are the reaction forces and R is
% the cholesky factorization of the free part of K.
%
% [u, Q, R, B] = MSOLVEQ(K, f, bc, Kold, R, s) approximates the solution u
% to the linear system of equations using CA. Kold is the stiffness matrix
% corresponding to the cholesky facotorization R. s orthonormal basis
% vectors are computed using CASBON
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
    [V, B] = CASBON(R, dK, Kff, ff, s);
    uf = V*(V'*ff);
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