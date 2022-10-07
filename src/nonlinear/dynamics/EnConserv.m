function [a, da] = EnConserv(Kt, M, fint, fext, h, IV, bc)
% Computes the dynamical response to the system with an energy conserving
% algorithm
% 
% Inputs:
%   Kt:     Stiffness matrix                (function with 2 inputs)
%
%   M:      Mass matrix                     [m x m] 
%
%   fint:   Internal force vector           (function with 2 inputs)
%
%   fext:   External force vector           [m x 1]
%
%   h:      Time step size
%
%   IV:     initial values                  [m x 2]
%
%   bc:     boundary conditions             [k x 2]
rtol = 1e-3;
utol = 1e-6;
[n, m] = size(fext);

alpha = (1:n);
if ~isempty(bc)
    alpha(bc(:, 1)) = [];
end

tend = m * h;
[updateBar, wbar] = makeMonitor(tend);

Keff = @(ap, dela) (2/h)^2*M + Kt(ap + dela, ap);
reff = @(ap, dela, dap, fmext) ...
    (2/h)^2*M*dela + 2*fint(ap + dela, ap) - 2*fmext - 4/h*M*dap;

a = zeros([n, m]);
a(:, 1) = IV(:, 1);

da = zeros([n, m]);
da(:, 1) = IV(:, 2);

for k = 2:m
    % Checking status of monitor (bar)
    if getappdata(wbar, 'canceling')
        break;
    end
    % Initiating variables
    fmext = (fext(:, k) + fext(:, k-1))/2;
    
    ap = a(:, k-1);     % Previous displacement
    dap = da(:, k-1);   % Previous velocity
    
    % Predictor
    dela = h*dap;
    
    % Inner newton scheme
    ninner = 0;
    rc = reff(ap, dela, dap, fmext);
    while ninner == 0 || (norm(rc(alpha)) > rtol && norm(Da) > utol)
        K = Keff(ap, dela);
        if ~isempty(bc)
            Da = solveq(K, -rc, bc);
        else
            Da = solveq(K, -rc);
        end
        % Updating
        dela = dela + Da;
        
        % Out of balance forces
        rc = reff(ap, dela, dap, fmext);
        
        ninner = ninner + 1;
        if ninner > 50
            delete(wbar)
            error("Too many iterations")
        end
        updateBar(k/m, h*k);
    end
    a(:, k) = ap + dela;
    da(:, k) = (2/h)*dela - dap;
end
delete(wbar)
end

function printOuter(k)
fprintf('%4i\n', k)
end

function printInner(rc)
fprintf('\t %4.2d\n', norm(rc))
end

function [updateBar, wbar] = makeMonitor(tend)
    wbar = waitbar(0, '1', 'Name', 'Energy conserving algorithm', ...
        'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', 1)');
    setappdata(wbar, 'canceling', 0);
    updateText = @(t) sprintf('Current time: %.2fms (tend: %ims)', ...
        t*1000, floor(tend*1000));
    updateBar = @(q, t) waitbar(q, wbar, updateText(t));
end