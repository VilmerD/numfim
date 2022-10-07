function [a, da, dda] = Newmark(Kt, d1, M, fint, fext, h, beta, gamma, ...
    IV, bc, predictor)
% Simulates the dynamical response of the system
%
% Inputs:
%   Kt:         Stiffness matrix                (function with 1 input)
%
%   M:          Mass matrix                     [m x m] 
%
%   fint:       Internal force vector           (function with 1 input)
%
%   fext:       External force vector           [m x 1]
%
%   h:          Time step size
%
%   beta:       Parameter which must fullfil 
%               0.5 < g < 1
%
%   gamma:      Parameter which must fullfil 
%               0.5*g < b < 1/2
%
%   IV:         initial values                  [m x 2]
%
%   bc:         boundary conditions             [k x 2]
%
%   predictor:  method used to predict next step
%               can be "Acc", "Vel", or "None"
%
% The method throws an error if the residual grows beyond 1e15.
if gamma < 0.5 || gamma > 1
    error("Bad choice of gamma, should be 0.5 < g < 1")
elseif beta < 0.5*gamma || beta > 1
    error("Bad choice of beta, should be 0.5*g < b < 1");
end

if nargin < 11
   predictor = 'None';
end

rtol = 1e-4;
utol = 1e-9;

[n, m] = size(fext);
s = [n, m];

alpha = (1:n);
if ~isempty(bc)
    alpha(bc(:, 1)) = [];
end

tend = m * h;
[updateBar, wbar] = makeMonitor(tend);

% Computing some constants used in the method
g = gamma;  b = beta;
c1 = 1/(b*h^2); % c2 = 1/(b*h);   c3 = (1 - 2*b)/2*b;       % Not needed
c4 = g/(b*h);   % c5 = (g - b)/b; c6 = h*(g - 2*b)/(2*b);

Keff = @(a) (c1 + d1*c4)*M + Kt(a);
reff = @(a, da, dda, fkext) M*dda + d1*M*da + fint(a) - fkext;

% Initiating solution vectors, and inserting boundary values
a = zeros(s);
a(:, 1) = IV(:, 1);
da = zeros(s);
da(:, 1) = IV(:, 2);
dda = zeros(s);
dda(:, 1) = -d1*da(:, 1) - solveq(M, (fint(a(:, 1)) - fext(:, 1)), bc);

for k = 2:m
    % Checking status of monitor (bar)
    if getappdata(wbar, 'canceling')
        break;
    end
    % Initiating variables
    fkext = fext(:, k);
    
    ak = a(:, k-1);
    dak = da(:, k-1);
    ddak = dda(:, k-1);
    
    % Predictor
    if strcmpi("Vel", predictor)
        [ak, dak, ddak] = constVel(ak, dak, ddak, beta, gamma, h);
    elseif strcmpi("Acc", predictor)
        [ak, dak, ddak] = constAcc(ak, dak, ddak, h);
    elseif ~strcmpi("None", predictor)
        error("Unrecognized predictor option, choices are: " ... 
            + "(constant) Vel or Acc.")
    end
    
    % Inner newton scheme
    ninner = 1;
    rc = reff(ak, dak, ddak, fkext);
    while ninner == 1 || (norm(rc(alpha)) > rtol && norm(Da) > utol)
        % Taking a step
        K = Keff(ak);
        Da = solveq(K, -rc, bc);
        
        % Updating quantities
        ak = ak + Da;
        dak = dak + c4*Da;
        ddak = ddak + c1*Da;
                
        % Out of balance forces
        rc = reff(ak, dak, ddak, fkext);
        
        ninner = ninner + 1;
        if norm(rc(alpha)) > 1e15
            delete(wbar)
            error("Too large residual")
        end
        updateBar(k/m, h*k);
    end
    
    % Accepting quantities
    a(:, k) = ak;
    da(:, k) = dak;
    dda(:, k) = ddak;
    
end
delete(wbar)
end

% Predictor of next step assuming constant velocity
function [an, dan, ddan] = constVel(ak, dak, ddak, b, g, h)
ddan = (1 - g)/g*ddak;
dan = dak;
an = ak + dak*h + ((1 - 2*b)*ddak + 2*b*ddan)*h^2/2;
end

% Predictor of next step assuming constant acceleration
function [an, dan, ddan] = constAcc(ak, dak, ddak, h)
ddan = ddak;
dan = dak + h * ddak;
an = ak + h*dak + h^2/2 * ddak;
end

% Creates the waitbar (monitor)
function [updateBar, wbar] = makeMonitor(tend)
    wbar = waitbar(0, '1', 'Name', 'Newmarks algorithm', ....
        'CreateCancelBtn', 'setappdata(gcbf, ''canceling'', 1)');
    setappdata(wbar, 'canceling', 0);
    updateText = @(t) sprintf('Current time: %.2fms (tend: %ims)', ...
        t*1000, floor(tend*1000));
    updateBar = @(q, t) waitbar(q, wbar, updateText(t));
end