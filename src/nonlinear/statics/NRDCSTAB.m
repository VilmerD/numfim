function [P, u] = NRDCSTAB(K, r, xp, nmax, ndof, options)
% Computes the displacement controlled response of the system
%
% Inputs
%   Kt:     Stiffness matrix        (function with 1 input)
%
%   r:      Residual force vector   (function with 2 inputs)
%
%   xp:     Prescribed displacement [nbc x 1]
%
%   nmax:   amount of steps
%   
%   ndof:   number of degrees of freedom
%
%   options:
%
rMax = 1e12;
if isfield(options, 'rtol')
    rtol = options.rtol;
else

rtol = 1e-7;
end

if isfield(options, 'u0')
    u0 = options.u0;
    n0 = options.n0;
else
    u0 = zeros(ndof, 1);
    n0 = 1;
end

if isfield(options, 'Verbose')
    verbosity = options.Verbose;
else
    verbosity = 1;
end

if isfield(options, 'N_INNER_MAX')
    N_INNER_MAX = options.N_INNER_MAX;
else
    N_INNER_MAX = 10;
end

% Reset warning so illConditionedMatrix warning can be caught
lastwarn('');
s = warning('off', 'MATLAB:illConditionedMatrix');

% The boundary condition used to correct the solution
correct = xp;
correct(:, 2) = 0;
resn = r(u0, 0);

% Free/Prescribed nodes
np = xp(:, 1);
nf = 1:ndof;
nf(np) = [];

% How much to load each step
dup_k = xp(:, 2)/nmax;            

% Initiating quantities
un = u0;
u = zeros(ndof, nmax);
P = zeros(ndof, nmax);
if verbosity
    printHeading();
end
for n = n0:nmax
    if verbosity
        fprintf('\n%12s%9i%9s', 'Taking Step', n, '');
    end
    % Computing how much should be displaced
    dup_n = n*dup_k - un(np);
    load = [xp(:, 1) dup_n];
    
    % Displacement
    if isfield(options, 'solver')
        [~, s] = options.solver(K(un), -resn, load, n);
    else
        s = solveq(K(un), -resn, load);
    end
    % Line search
    r0 = -s'*resn;
    r1 = -s'*r(un + s, 0);
    bapp = r0/(r0 - r1);
    b = max(0.3, min(bapp, 3));
    dun = s;
    dun(nf) = b*s(nf);
    
    un = un + dun;
    resn = r(un, 0);                % Residual foces
    r_free = norm(resn(nf));        % Norm of residual in free nodes
    
   % Check for warnings when assembling r
    checkResidualWarnings();
    
    fprintf('% 1.2e', r_free)
    
    % Iterating untill convergance
    N_INNER = 0;
    if verbosity
        fprintf('\n%12s%9i%9i%9s', 'Correcting', n, '', '');
    end
    while r_free > rtol
        % Computing new estimate, with zero displacement in prescribed nodes
        if isfield(options, 'solver')
            [~, s] = options.solver(K(un), -resn, correct, n);
        else
            s = solveq(K(un), -resn, correct);
        end
        % Line search
        r0 = -s'*resn;
        r1 = -s'*r(un + s, 0);
        bapp = r0/(r0 - r1);
        b = max(0.3, min(bapp, 3));
        dun = s;
        dun(nf) = b*s(nf);
        
        un = un + dun;
        resn = r(un, 0);
        r_free = norm(resn(nf));
        
        checkResidualWarnings();
        
        N_INNER = N_INNER + 1;
        % Update user
        if verbosity
            printAction('', n, N_INNER, r_free);
        end
        if N_INNER > N_INNER_MAX || r_free > rMax
            errorStruct.message = sprintf(...
                'Newton failed to converge within %i steps.', N_INNER_MAX);
            errorStruct.identifier = 'NR:ConverganceError';
            error(errorStruct)
        end
    end
    u(:, n) = un;
    P(:, n) = resn;
end
end

function printHeading()
fprintf('\n   Action    n_outer  n_inner     r    ');
end

function printAction(action, nouter, ninner, r)
fprintf('\n%12s%9i%9i% 1.2e', action, nouter, ninner, r)
end

function checkResidualWarnings()
    [warnmsg, ~] = lastwarn;
    if ~isempty(warnmsg)
        errorStruct.message = 'Problems with deformation gradient';
        errorStruct.identifier = 'NR:FE_Error';
        error(errorStruct);
    end
end