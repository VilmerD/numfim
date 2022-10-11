function [P, L, d, V_INTR] = CAeigs(A, B, n, Rold, Bold, Pold, s, options)
% CAeigs(A, B, n, Rold, Bold, Pold, s, options) finds the n eigenparis of the
% generalized eigenvalue problem correspodning to the n eigenvalues with
% the largest absolute value.
% 
% R is the cholesky factorization of Bold, an older version of B 
% Pold are the old eigenvectors and s is the number of basis vectors to use.
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

P = zeros(size(A, 1), n);
L = zeros(n, n);
d = -1*ones(n, 1);
V_INTR = cell(n, 4);
for k = 1:n
    % Choose if orthogonalization is to be used and which vectors to
    % orthogonalize with respect to
    switch upper(options.orthotype)
        case 'CURRENT'
            VO = P(:, 1:(k-1));
        case 'OLD'
            VO = options.orthovecs(:, 1:(k - 1));
        otherwise
            VO = [];
    end
    
    % Generate basis vectors
    [Vk, Uk, Tk] = CAEEON(A, B, Rold, Bold, Pold(:, k), s, VO);
    
    % Compute reduced model
    A_RED = Vk'*A*Vk;
    B_RED = Vk'*B*Vk;
    
    % Solve reduced problem
    [P_RED, L_RED] = eigs(A_RED, B_RED, 1, 'largestabs', ...
        'IsCholesky', false);
    
    % Normalize the eigenvectors with respect to A
    P_RED = P_RED/sqrt(P_RED'*B_RED*P_RED);                    
    P_FULL = Vk*P_RED;
    
    % Insert solution
    P(:, k) = P_FULL;                             
    L(k, k) = L_RED;
    V_INTR(k, :) = {Uk, Tk, Vk, P_RED};
    
    % Compute the residual of the eigenproblem
    d(k) = norm(A*P_FULL - L_RED*B*P_FULL)/norm(A*P_FULL);
end

end