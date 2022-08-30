function dw = sCAERS(B, psi, w, K, M, Kold, R, psi0, edof, np, K0, M0)
% sCAERS computes the senisitivy of the eigenvalues when Rayleigh-shifted
% CA is used to reduce the eigenproblem

ndof = size(K, 1);
nf = (1:ndof)';
nf(np) = [];
Kff = K(nf, nf);
Mff = M(nf, nf);

dKff = Kff - Kold(nf, nf);

% Allocating memory
ne = numel(w);
Zs = cell(ne, 1);
QUs = cell(ne, 1);
uns = cell(ne, 1);
muks = cell(ne, 1);
rmus = cell(ne, 1);
yhats = cell(ne, 1);

% Computing adjoint vectors
nelm = size(edof, 1);
for k = 1:ne
    yk = B{k, 2};
    Ufk = B{k, 1}{1};
    Vfk = B{k, 1}{2};
    s = size(yk, 1);
    
    psif = psi(nf, k);
    muk = sqrt(psif'*Mff*psif);
    dff = (Kff*psif - w(k)^2*Mff*psif)/muk^2;
    
    % Computing the LAST vector first
    qs = -2*dff*yk(s);
    
    us = Ufk(:, s);
    qsus = qs'*us;
    
    muf = Mff*Ufk;
    mus = muf(:, s);
    usn = sqrt(us'*mus);            
    Zk = zeros(ndof, s);
    Zk(nf, s) = R\(R'\((usn^2*qs - qsus*mus)))/usn^3;
    
    QU = zeros(s, 1);
    QU(s) = qsus;
    un = zeros(s, 1);               % M-norm of u-vectors
    un(s) = usn;                    
    
    mvf = Mff*Vfk;
    mu = zeros(s, 1);               % Shift constants
    yhat = zeros(s, 1);
    % Stepping backward
    for i = (s-1):-1:1
        vi = Vfk(:, i);
        mvi = mvf(:, i);
        muip1 = (vi'*Kff*vi)/(vi'*mvi);
        dfmuip1 = (Kff*vi - muip1*mvi);
        zip1 = Zk(nf, i + 1);
        yhatip1 = vi'*Mff*zip1;
        qi = -2*yk(i)*dff - dKff*zip1 + muip1*Mff*zip1 + 2*dfmuip1*yhatip1;
        
        ui = Ufk(:, i);
        mui = muf(:, i);
        qiui = qi'*ui;
        uin = sqrt(ui'*mui);
        Zk(nf, i) = R\(R'\((uin^2*qi - qiui*mui)))/uin^3;
        
        QU(i) = qiui;
        un(i) = uin;
        mu(i + 1) = muip1;
        yhat(i + 1) = yhatip1;
    end
    Zs{k} = Zk;
    QUs{k} = QU;
    uns{k} = un;
    muks{k} = muk;
    rmus{k} = mu;
    yhats{k} = yhat;
end

Us = cell(ne, 1);
Vs = cell(ne, 1);
for i = 1:ne
    Bi = B{i, 1};
    si = size(Bi{1}, 2);
    Usi = zeros(ndof, si);
    Vsi = Usi;
    Usi(nf, :) = Bi{1};
    Vsi(nf, :) = Bi{2};
    
    Us{i} = Usi;
    Vs{i} = Vsi;
end

% Computing sensitivity on the element level
dl = zeros(ne, nelm);
for elm = 1:nelm
    % Preparing matrices
    dofs = edof(elm, 2:end);
    
    ke = K0{elm};
    me = M0{elm};
    
    % Computing sensitivity of each eigenvalue
    for k = 1:ne
        Zelm = Zs{k}(dofs, :);
        Uelm = Us{k}(dofs, :);
        Velm = Vs{k}(dofs, :);
        sk = size(Zelm, 2);
        QU = QUs{k};
        un = uns{k};
        muk = muks{k};
        mu = rmus{k};
        yhatk = yhats{k};
        
        psie = psi(dofs, k);
        psi0e = psi0(dofs, k);
        dlelm = psie'*(ke - w(k)^2*me)*psie/muk^2 - (Zelm(:, 1))'*me*psi0e;
        
        u1elm = Uelm(:, 1);
        dlelm = dlelm + 1/2*QU(1)*(u1elm'*me*u1elm)/un(1)^3;
        
        for i = 2:sk
            zielm = Zelm(:, i);
            vim1 = Velm(:, i-1);
            smui = (ke - mu(i)*me);
            dlelm = dlelm + zielm'*smui*vim1 - yhatk(i)*(vim1'*smui*vim1);
            
            uielm = Uelm(:, i);
            dlelm = dlelm + 1/2*QU(i)*(uielm'*me*uielm)/un(i)^3;
        end
        dl(k, elm) = dlelm;
    end
end
end