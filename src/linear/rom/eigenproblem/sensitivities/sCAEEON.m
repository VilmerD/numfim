function dldz = sCAEEON(B, psi, L, K, M, Kold, R, psi0, edof, np, K0, M0)
% sCAEEON computes the sensitivity of the eigenvalues when
% eigenmode orthogonalized CA is used to reduce the eigenproblem
ndof = size(K, 1);
nf = (1:ndof)';
nf(np) = [];
Kff = K(nf, nf);
Mff = M(nf, nf);

dKff = Kff - Kold(nf, nf);

% Preallocating memory
nelm = size(edof, 1);
ne = numel(L);
Zs = cell(ne, 1);
Us = cell(ne, 1);
Ts = cell(ne, 1);

QU_UNs = cell(ne, 1);
WPs = cell(ne, 1);

mpsi = M*psi;
% Computing adjoint vectors
for k = 1:ne
    % Extracting data
    Bk = B{k, 1};
    yk = B{k, 2};
    sk = size(yk, 1);
    Uk = zeros(ndof, sk);
    Tk = zeros(ndof, sk);
    Uk(nf, :) = Bk{1};
    Tk(nf, :) = Bk{2};
    
    psikf = psi(nf, k);
    dff = (Kff*psikf - L(k)*Mff*psikf);
    
    % Computing the LAST vector first
    ws = -2*dff*yk(sk);
    
    % last q-adjoint
    qs = ws;
    WP = zeros(sk, k-1);
    for j = 1:(k - 1)
        wspj = ws'*psi(nf, j);
        WP(sk, j) = wspj;
        qs = qs - mpsi(nf, j)*wspj;
    end
    
    mu = M*Uk;
    us = Uk(nf, sk);
    mus = mu(nf, sk);
    usn = sqrt(us'*mus);
    qsus_usn = qs'*us/usn^3;
    
    Zk = zeros(ndof, sk);
    Zk(nf, sk) = R\(R'\((qs/usn - qsus_usn*mus)));
    
    QU_UN = zeros(sk, 1);
    QU_UN(sk) = qsus_usn;
    
    % Stepping backward
    for i = (sk-1):-1:1
        wi = -2*yk(i)*dff;
        qi = wi - dKff*Zk(nf, i+1);
        
        for j = 1:(k - 1)
            wipj = wi'*psi(nf, j);
            qi = qi - mpsi(nf, j)*wipj;
            WP(i, j) = wipj;
        end
        
        ui = Uk(nf, i);
        uin = sqrt(ui'*mu(nf, i));
        qiui_uin = qi'*ui/uin^3;
        Zk(nf, i) = R\(R'\((qi/uin - (mu(nf, i))*(qiui_uin))));
        
        QU_UN(i) = qiui_uin;
    end
    
    Zs{k} = Zk;
    Us{k} = Uk;
    Ts{k} = Tk;
    QU_UNs{k} = QU_UN;
    WPs{k} = WP;
end

% Index matrix to extract element data
dldz = zeros(ne, nelm);
% Computing sensitivity on the element level
for elm = 1:nelm
    % Preparing matrices
    dofs = edof(elm, 2:end);
    psi0e = psi0(dofs, :);
    psie = psi(dofs, :);
    
    ke = K0{elm};
    me = M0{elm};
    
    % Computing sensitivity
    for k = 1:ne
        Zelm = Zs{k}(dofs, :);
        Uelm = Us{k}(dofs, :);
        Telm = Ts{k}(dofs, :);
        QU_UN = QU_UNs{k};
        WP = WPs{k};
        
        psiek = psie(:, k);
        psi0ek = psi0e(:, k);
        
        dlelm = psiek'*(ke - L(k)*me)*psiek;
        
        zterm = - (Zelm(:, 1))'*me*psi0ek;
        
        qterm = 1/2*QU_UN(1)*(Uelm(:, 1)'*me*Uelm(:, 1));
        for i = 2:sk
            zterm = zterm + (Zelm(:, i))'*ke*Telm(:, i-1);
            
            qterm = qterm + 1/2*QU_UN(i)*(Uelm(:, i)'*me*Uelm(:, i));
        end
        
        wterm = 0;
        if k > 1
            wterm = sum(WP.*(Telm'*me*psie(:, 1:(k - 1))), 'all');
        end
        
        dldz(k, elm) = dlelm + zterm + qterm + wterm;
    end
    
end
end