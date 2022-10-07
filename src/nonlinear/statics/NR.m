function [P, u, N] = NR(Kt, f, p, nmax, bc, u0)
% Computes the load controlled response of the system
%
% Inputs
%   Kt:     Stiffness matrix        (function with 1 input)
%
%   f:      Internal force vector   (function with 1 input)
%
%   p:      Final force             [m x 1]
%
%   nmax:   amount of steps
%   
%   bc:     boundary conditions     [k x 2]
%
%   u0:     initial displacement    [m x 1]
%           optional, default is 0
rtol = 1e-6;
dP = p/nmax;

if nargin < 6
    u0 = zeros(length(dP), 1);
end
m = length(dP);

alpha = (1:m);
if ~isempty(bc)
    alpha(bc(:, 1)) = [];
end

[updateBar, wbar] = makeMonitor(nmax);

P = zeros(length(dP), nmax + 1);
u = [u0, zeros(length(u0), nmax)];
NUMBER_OF_ITERATIONS = 0;
for n = 2:nmax + 1
    % Checking status of monitor (bar)
    if getappdata(wbar, 'canceling')
        break;
    end
    Pn = P(:, n-1) + dP;
    P(:, n) = Pn;

    un = u(:, n - 1);
    rc = f(un) - Pn;
    
    nit_inner = 0;
    while nit_inner == 0 || norm(rc(alpha)) > rtol
        K = Kt(un);

        [du, ~] = solveq(K, -rc, bc);
        un = un + du;

        rc = f(un) - Pn;

        nit_inner = nit_inner + 1;
        if nit_inner > 40
            delete(wbar)
            error("Too large residual")
        end
    end
    u(:, n) = un;
    updateBar(n/nmax, n);
    NUMBER_OF_ITERATIONS = NUMBER_OF_ITERATIONS + nit_inner;
end
N = NUMBER_OF_ITERATIONS;
delete(wbar)
end


function [updateBar, wbar] = makeMonitor(nmax)
    wbar = waitbar(0, '1', 'Name', 'Load controlled scheme', ...
        'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', 1)');
    setappdata(wbar, 'canceling', 0);
    updateText = @(n) sprintf(['Current step: %i (nmax: %i)'], n, nmax);
    updateBar = @(q, n) waitbar(q, wbar, updateText(n));
end