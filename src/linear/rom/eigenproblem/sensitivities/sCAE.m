function dldz = sCAE(B, psi, L, K, M, Kold, R, psi0, edof, np, K0, M0)
% sCAE computes the sensitivity of the eigenvalues when CA is used to
% reduce the eigenproblem
ndof = size(K, 1);
nf = (1:ndof)';
nf(np) = [];
Kff = K(nf, nf);
Mff = M(nf, nf);

dKff = Kff - Kold(nf, nf);

% Allocating memory
ne = numel(L);
Zs = cell(ne, 1);
QUs = cell(ne, 1);
uns = cell(ne, 1);

% Computing adjoint vectors
for k = 1:ne
    yk = B{k, 2};
    Ufk = B{k, 1}{1};
    s = size(yk, 1);
    
    psif = psi(nf, k);
    dff = (Kff*psif - L(k)*Mff*psif);
    
    % Computing the LAST vector first
    qs = -2*dff*yk(s);
    
    us = Ufk(:, s);
    qsus = qs'*us;
    
    muf = Mff*Ufk;
    mus = muf(:, s);
    usn = sqrt(us'*mus);
    Zk = zeros(ndof, s);
    Zk(nf, s) = R\(R'\((usn^2*qs - qsus*mus)))/usn^3;
    
    Qk = zeros(ndof, s);
    Qk(nf, s) = qs;
    
    QU = zeros(s, 1);
    QU(s) = qsus;
    un = zeros(s, 1);
    un(s) = usn;
    
    % Stepping backward
    for i = (s-1):-1:1
        qi = -2*yk(i)*dff - dKff*Zk(nf, i+1);
        
        ui = Ufk(:, i);
        mui = muf(:, i);
        qiui = qi'*ui;
        uin = sqrt(ui'*mui);
        Zk(nf, i) = R\(R'\((uin^2*qi - qiui*mui)))/uin^3;
        
        Qk(nf, i) = qi;
        QU(i) = qiui;
        un(i) = uin;
    end
    Zs{k} = Zk;
    QUs{k} = QU;
    uns{k} = un;
end

% Create 'full' basis vectors
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
nelm = size(edof, 1);
edof_red = edof(:, 2:end);
dldz = zeros(ne, nelm);
for elm = 1:nelm
    % Preparing matrices
    dofs = edof_red(elm, :);
    
    ke = K0{elm};
    me = M0{elm};
    psie = psi(dofs, :);
    psi0e = psi0(dofs, :);
    
    % Computing sensitivity of each eigenvalue
    for k = 1:ne
        Zelm = Zs{k}(dofs, :);
        Uelm = Us{k}(dofs, :);
        Velm = Vs{k}(dofs, :);
        sk = size(Zelm, 2);
        QU = QUs{k};
        un = uns{k};
        
        psiek = psie(:, k);
        psi0ek = psi0e(:, k);
        dlelm = psiek'*(ke - L(k)*me)*psiek - (Zelm(:, 1))'*me*psi0ek;
        
        u1elm = Uelm(:, 1);
        dlelm = dlelm + 1/2*QU(1)*(u1elm'*me*u1elm)/un(1)^3;
        
        for i = 2:sk
            dlelm = dlelm + (Zelm(:, i))'*ke*Velm(:, i-1);
            
            uielm = Uelm(:, i);
            dlelm = dlelm + 1/2*QU(i)*(uielm'*me*uielm)/un(i)^3;
        end
        dldz(k, elm) = dlelm;
    end
end
end