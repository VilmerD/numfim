function [P, u, N] = NRLC(Kt, r, p, nmax, bc, options)
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

% Setup some data for the solver
if isfield(options, 'rtol')
    rtol = options.rtol;
else
    rtol = 1e-4;
end
m = length(p);

% Number of maximum inner iterations, if number of iteration surpasses this
% number an error is thrown
if isfield(options, 'N_INNER_MAX')
    N_INNER_MAX = options.N_INNER_MAX;
else
    N_INNER_MAX = 50;
end

% Reset warning so illConditionedMatrix warning can be caught
lastwarn('');
s = warning('off', 'MATLAB:illConditionedMatrix');

% Initial guess
if isfield(options, 'u0')
    u0 = options.u0;
    n0 = options.n0;
else
    u0 = zeros(m, 1);
    n0 = 1;
end
nsteps = (nmax - (n0 - 1));
dP = p/nsteps;

% Setup boundary condition and free noes
nf = (1:m);
np = bc(:, 1);
if ~isempty(bc)
    nf(np) = [];
end

% Initialize solution vectors
P = zeros(m, nmax + 1);
u = zeros(m, nmax + 1);
u(:, n0) = u0;
NUMBER_OF_ITERATIONS = 0;
for n = n0:nmax
    % Take a load step
    Pn = P(:, n) + dP;
    P(:, n + 1) = Pn;
    
    % Initialize displacement vector and residual
    un = u(:, n);
    rn = r(un, Pn);
    r_free = norm(rn(nf));
    
    % Check for warnings when assembling r
    [warnmsg, ~] = lastwarn;
    if ~isempty(warnmsg)
        errorStruct.message = 'Problems with deformation gradient';
        errorStruct.identifier = 'NR:FE_Error';
        error(errorStruct);
    end
    
    n_inner = 0;
    % Correction step
    while r_free > rtol
        K = Kt(un);

        % Solve for correction and update
        if isfield(options, 'solver')
            % Generalize this so the state of the iteration is sent as a
            % parameter
            [~, du] = options.solver(K, -rn, bc, n);
        else
            [~, du] = solveq(K, -rn, bc);
        end
        
        % Update quantities
        un = un + du;
        rn = r(un, Pn);
        r_free = norm(rn(nf));

        n_inner = n_inner + 1;
        fprintf('\n(NR) r: %1.2e', r_free);
        if n_inner > N_INNER_MAX
            errorStruct.message = sprintf(...
                'Newton failed to converge within %i steps.', N_INNER_MAX);
            errorStruct.identifier = 'NR:ConverganceError';
            error(errorStruct)
        end
    end
    
    % Accept quantities
    u(:, n + 1) = un;
    NUMBER_OF_ITERATIONS = NUMBER_OF_ITERATIONS + n_inner + 1;
    
end
N = NUMBER_OF_ITERATIONS;

% Reset warning
warning(s);
end