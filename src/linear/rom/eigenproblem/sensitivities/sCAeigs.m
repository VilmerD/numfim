function DL_Dzf = sCAeigs(B, psi, L, K, M, Kold, R, psi0, ...
    edof, np, DK_Dzf, DM_Dzf, options)
% Generic method, i don't ork explain now

switch options.orthotype
    case 'none'
        DL_Dzf = sCAE(B, psi, L, K, M, Kold, R, psi0, edof, np, ...
            DK_Dzf, DM_Dzf);
    case 'current'
        DL_Dzf = sCAEEON(B, psi, L, K, M, Kold, R, psi0, edof, np, ...
            DK_Dzf, DM_Dzf);
    otherwise
        vectors = options.orthovecs;
        DL_Dzf = sCAEEONmod(B, psi, L, K, M, Kold, R, psi0, vectors, ...
            edof, np, DK_Dzf, DM_Dzf);
end

end