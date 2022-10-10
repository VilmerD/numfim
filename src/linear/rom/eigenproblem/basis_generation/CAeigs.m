function [V, D, deltas, B] = CAeigs(K, M, n, R, dK, psi0, s, options)
% CAeigs(K, M, n, R, dK, psi0, s, options) finds the n eigenparis of the
% generalized eigenvalue problem correspodning to the n eigenvalues with
% the smallest absolute value.
% 
% R is the cholesky factorization, dK is the change in stiffness matrix,
% psi0 are the old eigenvectors and s is the number of basis vectors to use.
% 
% Finally options is a struct containing the orthotype and orthovecs
% 
% If orthotype is CURRENT CA uses the current approximation of the
% eigenvectors for orthogonalization.
% 
% If orthotype is OLD CA uses an old approixmation of the eigenvectors
% for orthogonalization, which must be contained in orthovecs
%
% If orthotype is NONE or otherwise, CA does not do basis
% orthogonalization.

V = zeros(size(K, 1), n);
D = zeros(n, n);
deltas = -1*ones(n, 1);
B = cell(n, 2);
for k = 1:n
    % Choose if orthogonalization is to be used and which vectors to
    % orthogonalize with respect to
    switch upper(options.orthotype)
        case 'CURRENT'
            vectors = V(:, 1:(k-1));
        case 'OLD'
            vectors = options.orthovecs(:, 1:(k - 1));
        otherwise
            vectors = [];
    end
    
    % Generate basis vectors
    [Vk, Bk] = CAEEON(R, dK, K, M, psi0(:, k), vectors, s);    
    % Compute reduced model
    KR = Vk'*K*Vk;
    MR = Vk'*M*Vk;
    
    % Solve reduced model
    [yk, dk] = eigs(KR, MR, 1, 'smallestabs', ...
        'IsCholesky', false);
    
    % Normalize the eigenvectors with respect to the mass matrix
    yk = yk/sqrt(yk'*MR*yk);                    
    psik = Vk*yk;
    
    % Insert solution
    V(:, k) = psik;                             
    D(k, k) = dk;
    B{k, 1} = Bk;
    B{k, 2} = yk;
    
    % Compute the residual of the eigenproblem
    deltas(k) = norm(K*psik - dk*M*psik)/norm(K*psik);
end

end