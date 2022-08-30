function [u, P, ef, es, flagout] = NRDCoptstep(K, r, sfun, dz, zold, bc, ...
    uold, RESTARTS_MAX, NSTEPS, options)
fprintf('%s\n', repmat('-', 1, 33));

% Initialize flags and counters
PROBLEM_SOLVED = 0;
NRESTARTS = 0;
START_AT_ZERO = 0;

% Initialize stiffness matrix etc
znew = zold + dz;
Kz = @(ef, es) K(ef, es, znew);
rz = @(ef, es, fext) r(ef, es, fext, znew);

% Restarts NRDC at steps u0
while ~PROBLEM_SOLVED
    
    % Get starting point
    if ~START_AT_ZERO && NRESTARTS <= RESTARTS_MAX
        u0n = uold(:, end - NRESTARTS);
        fprintf('Start number %i\n', NRESTARTS);
    elseif START_AT_ZERO
        fprintf('Starting from 0\n');
        u0n = zeros(size(uold, 1), 1);
        options.solver.forceFactorization = 1;
    else
        errorStruct.message = 'Iteration failed';
        error(errorStruct);
    end
    
    % Attempt a solution at u0n
    [u, P, ef, es, NR_FLAG, NSTEP, N_INNER, relres] = ...
        NRDC(Kz, rz, sfun, bc, u0n, NSTEPS, options);
    
    % Check results
    if NR_FLAG == 0
        % Solution found! Get out!
        flagout = 0;
        PROBLEM_SOLVED = true;
        
    elseif NR_FLAG == 1 && NRESTARTS < RESTARTS_MAX
        % Restart with other initial guess
        NRESTARTS = NRESTARTS + 1;
        
    elseif NR_FLAG == 1 && NRESTARTS == RESTARTS_MAX
        % If dz is large it's feasible to divide into smaller parts,
        % otherwise, of if division fails start at zero
        if norm(dz) > 1e-1
            % Try dividing dz into smaller parts (stepping in z-dimension?)
            nzsteps = 8;
            alph = 1/nzsteps;
            fprintf('Dividing dz into %i steps \n', nzsteps);
            
            dzalph = dz*alph;
            uoldalph = uold(:, end);
            for i = 1:nzsteps
                fprintf('\t Step %i \n', i);
                zi = zold + i*dzalph;
                Kzi = @(ef, es) K(ef, es, zi);
                rzi = @(ef, es, fext) r(ef, es, fext, zi);
                [u, P, ef, es, NR_FLAG] = NRDC(Kzi, rzi, sfun, bc, ...
                    uoldalph, NSTEPS, options);
                uoldalph = u(:, end);

                % If subproblem wasn't solved, break
                if NR_FLAG == 1
                    break;
                end
            end
        end
        
        % If dz is too small or division failed, start at zero
        if norm(dz) <= 1e-1 || NR_FLAG == 1
            START_AT_ZERO = true;
        else
            flagout = 0;
            PROBLEM_SOLVED = true;
        end
        
    elseif START_AT_ZERO
        % Started from zero, and coulnd't solve problem
        flagout = 2;
        fprintf('Failed when starting at zero \n');
        break;
    end
    
end

% If solution is still not found there is a problem...
if ~PROBLEM_SOLVED
    errorStruct.message = 'Couldnt solve problem';
    error(errorStruct);
end

fprintf('Solution found!\n')
fprintf('%s\n', repmat('-', 1, 33));
end