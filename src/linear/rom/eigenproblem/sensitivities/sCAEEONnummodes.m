function dw = sCAEEONnummodes(B, psi, w, K, M, Kold, R, psi0, edof, np, K0, M0, dpsi)
% sCAEEON computes the sensitivity of the eigenvalues when CA is used to
% approximate the eigenproblem
ndof = size(K, 1);
nf = (1:ndof)';
nf(np) = [];
Kff = K(nf, nf);
Mff = M(nf, nf);

dKff = Kff - Kold(nf, nf);

% Computing adjoint vectors
nelm = size(edof, 1);
n = numel(w);
dw = zeros(nelm, n);

mpsi = zeros(ndof, n);
for k = 1:n
    mpsi(nf, k) = Mff*psi(nf, k);
end

for k = 1:n
    Bk = B{k, 1};
    yk = B{k, 2};
    s = size(yk, 1);
    Uk = zeros(ndof, s);
    Tk = zeros(ndof, s);
    Uk(nf, :) = Bk{1};
    Tk(nf, :) = Bk{2};
    
    psikf = psi(nf, k);
    muk = sqrt(psikf'*Mff*psikf);
    dff = (Kff*psikf - w(k)^2*Mff*psikf)/muk^2;
    
    % Computing the LAST vector first
    ws = -2*dff*yk(s);
    W = zeros(ndof, s);
    W(nf, s) = ws;
    
    % last q-adjoint
    qs = ws;
    WP = zeros(s, k-1);
    for j = 1:(k - 1)
        wspj = ws'*psi(nf, j);
        WP(s, j) = wspj;
        qs = qs - mpsi(nf, j)*wspj;
    end
    
    us = Uk(nf, s);
    qsus = qs'*us;
    
    % Compute M*u in advance to save some time
    mu = M*Uk;
    mus = mu(nf, s);
    usn = sqrt(us'*mus);
    Zk = zeros(ndof, s);
    Zk(nf, s) = R\(R'\((usn^2*qs - qsus*mus)))/usn^3;
    
    QU = zeros(s, 1);
    QU(s) = qsus;
    un = zeros(s, 1);
    un(s) = usn;
    
    % Stepping backward
    for i = (s-1):-1:1
        wi = -2*yk(i)*dff;
        qi = wi - dKff*Zk(nf, i+1);
        
        for j = 1:(k - 1)
            wipj = wi'*psi(nf, j);
            qi = qi - mpsi(nf, j)*wipj;
            WP(i, j) = wipj;
        end
        
        ui = Uk(nf, i);
        qiui = qi'*ui;
        uin = sqrt(ui'*mu(nf, i));
        Zk(nf, i) = R\(R'\((uin^2*qi - (mu(nf, i))*(qiui))))/uin^3;
        
        QU(i) = qiui;
        un(i) = uin;
        W(nf, i) = wi;
    end
    
    % Computing sensitivity on the element level
    dl = zeros(nelm, 1);
    for elm = 1:nelm
        % Preparing matrices
        dofs = edof(elm, 2:end);
        
        Zelm = Zk(dofs, :);
        Uelm = Uk(dofs, :);
        Telm = Tk(dofs, :);
        
        ke = K0{elm};
        me = M0{elm};
        
        % Computing sensitivity
        psie = psi(dofs, k);
        psi0e = psi0(dofs, k);
        dlelm = psie'*(ke - w(k)^2*me)*psie/muk^2 - (Zelm(:, 1))'*me*psi0e;
        
        for i = 2:s
            dlelm = dlelm + (Zelm(:, i))'*ke*Telm(:, i-1);
        end
        
        for i = 1:s
            uielm = Uelm(:, i);
            dlelm = dlelm + 1/2*QU(i)*(uielm'*me*uielm)/un(i)^3;
        end
        
        for i = 1:s
            for j = 1:(k - 1)
                psije = psi(dofs, j);
                dlelm = dlelm + WP(i, j)*(Telm(:, i)'*me*psije);
            end
        end
        
        % Adding numerically evaluated mode derivatives
        for i = 1:s
            Mti = M*Tk(:, i);
            wi = W(:, i);
            for j = 1:(k - 1)
                psij = psi(:, j);
                dpsijelm = dpsi(:, j, elm);
                dlelm = dlelm + ((Mti'*psij)*(wi'*dpsijelm) + WP(i, j)*(Mti'*dpsijelm));
            end
        end
        
        dl(elm) = dlelm;
    end
    
    dw(:, k) = dl/(2*w(k));
end
end