function [u, P, efn, esn, varargout] = NRDC(Kfun, rfun, sfun, bc, u0, NMAX, ...
    linear_solver)
% NRDC solves the nonlinear problem rfun = 0 with initial guess u0 using a
% newton-raphson scheme

% Computes the displacement controlled response of the system
ABS_MAX = 1e12;
N_INNER_MAX = 10;
RTOL = 1e-9;

% Free/Prescribed nodes
np = bc(:, 1);
nf = 1:size(u0, 1);
nf(np) = [];

% How much to load each step
[N_LOAD_STEPS, dustep] = computeNumberSteps(bc, u0, NMAX);
load_increment = [np dustep];
load_correct = [np zeros(size(np))];

% Initiating quantities
un = u0;
u = zeros(size(u0, 1), NMAX);
P = zeros(size(u0, 1), NMAX);

N_INNER_TOT = 0;
[efn, esn] = sfun(u0);
resn = rfun(efn, esn);
for n = (NMAX - N_LOAD_STEPS + 1):NMAX
    % Displacement increment
    Kn = Kfun(efn, esn);
    dun = linear_solver(Kn, -resn, load_increment, n);
    N_INNER = 1;
    
    % Updating displacement field and computing reaction forces
    un = un + dun;
    [efn, esn] = sfun(un);
    resn = rfun(efn, esn);
    rfree = norm(resn(nf));
    rtot = norm(resn);
    
    % Correcting until convergance
    while rfree/rtot > RTOL
        % Computing new estimate
        Kn = Kfun(efn, esn);
        dun = linear_solver(Kn, -resn, load_correct, n);
        N_INNER = N_INNER + 1;
        
        % Updating displacement field and computing reaction forces
        un = un + dun;
        [efn, esn] = sfun(un);
        resn = rfun(efn, esn);
        rfree = norm(resn(nf));
        rtot = norm(resn);
        
        % Checking if the iteration should be terminated
        if N_INNER > N_INNER_MAX || rfree > ABS_MAX || ~isreal(resn)
            if norm(u0) == 0
                % Started from zero and couldn't converge, bad
                flag = 2;
            else
                % Started from some initial guess, not as bad
                flag = 1;
            end
            
            % Packaging data
            N_INNER_TOT = N_INNER_TOT + N_INNER;
            argout = {flag, N_INNER_TOT, N_LOAD_STEPS, rfree/rtot};
            nargout_extra = nargout - 4;
            varargout = argout(1:nargout_extra);
            return
        end
    end
    % Inserting converged quantities
    u(:, n) = un;
    P(:, n) = resn;
    N_INNER_TOT = N_INNER_TOT + N_INNER;
end

% Packaging data
flag = 0;
argout = {flag, N_INNER_TOT, N_LOAD_STEPS, rfree/rtot};
nargout_extra = nargout - 4;
varargout = argout(1:nargout_extra);
end

function [nstep, step_new] = computeNumberSteps(bc, u0, nmax)
% The idea is to compute the number of steps it would take to reach bc from
% u0 if the same stepsize as if starting from 0 (default).
% This is done by first finding the total stepsize, and dividing it by the
% defualt stepsize. It is assumed that the load is applied linearly, ie by
% a set stepsize. Thus any nonzero element in bc can be chosen
np = bc(:, 1);          % Prescribed nodes

% Total load in:
up = bc(:, 2);          % prescribed nodes
upnz = up(up ~= 0);     % nonzero prescribed nodes
upnz1 = upnz(1);        % first nonzero prescribed node

% Initial disp in:
u0p = u0(np);           % prescribed nodes
u0pnz = u0p(u0p ~= 0);  % nonzero prescribed nodes

% first nonzero prescribed node
if isempty(u0pnz)
    u0pnz1 = 0;
else
    u0pnz1 = u0pnz(1);
end

% Additional displacement needed
du_add1 = upnz1 - u0pnz1;

% Default step size in first nonzero prescribed node
step_default1 = upnz1/nmax;

% Default stepsize compared to additional displacement
quo = round(du_add1/step_default1, 0);
nstep = max(1, quo);            % Ensure at least one step is taken
du_add = up - u0p;
step_new = du_add/nstep;
end