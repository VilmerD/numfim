function [x, b, varargout] = msolveq(A, b, bc, Rold, Aold, s)
% [x, b, R] = MSOLVEQ(A, b, bc) solves the linear system of equations 
% A*x = b with boundary conditions bc. Q are the reaction forces and R is
% the cholesky factorization of the free part of K.
%
% [x, b, R] = MSOLVEQ(A, b, bc, R) solves the linear system of equations 
% where R is the cholesky factorization of A.
%
% [x, b, U, T, R, V] = MSOLVEQ(A, b, bc, Rold, Aold, s) approximates the 
% solution u to the linear system of equations using CA. Kold is the 
% stiffness matrix corresponding to the cholesky facotorization R. s 
% orthonormal basis vectors are computed using CASBON
%

% If bc is empty just solve directly
if isempty(bc)
    Rcurr = chol(A);
    x = Rcurr\(Rcurr'\b);
    return
end

% Split up the system into free and prescribed nodes
ndof = size(A, 1);
nf = (1:ndof)';
np = bc(:, 1);
nf(np) = [];

Aff = A(nf, nf);

xp = bc(:, 2);
bf = b(nf) - A(nf, np)*xp;

% Solving free system
if nargin <= 4
    % Exactly
    if nargin < 4
        Rcurr = chol(Aff);
    else
        Rcurr = Rold;
    end
    xf = Rcurr\(Rcurr'\bf);
    varargout = {Rcurr};
else
    % Using CA
    [V, U, T, R] = CASBON(Aff, bf, Rold, Aold(nf, nf), s);
    xf = V*(V'*bf);
    varargout = {U, T, R, V};
end

% Reassembling the solution
x = zeros(ndof, 1);
x(np) = xp;
x(nf) = xf;

% Prescribed forces and reaction forces
b = zeros(ndof, 1);
b(np) = A(np, nf)*xf + A(np, np)*xp;
b(nf) = bf;
end