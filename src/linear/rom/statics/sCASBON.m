function dc = sCASBON(B, f, K, Kold, R, edof, np, K0)
% sCASBON computes the sensitivities of the compliance when CA is
% considered
nelm = size(edof, 1);
ndof = size(K, 1);
nf = (1:ndof)';
nf(np) = [];

Kff = K(nf, nf);
dKff = Kff - Kold(nf, nf);
ff = f(nf);

Uf = B{1};
Tf = B{2};
Rf = B{3};
Vf = B{4};

s = size(Uf, 2);
U = zeros(ndof, s);
T = zeros(ndof, s);
Rb = zeros(ndof, s);
V = zeros(ndof, s);
U(nf, :) = Uf;
T(nf, :) = Tf;
Rb(nf, :) = Rf;
V(nf, :) = Vf;

% Computing adjoint vectors
ps = -2*ff*(Vf(:, s)'*ff);
rs = Rb(nf, s);
nrs = sqrt(rs'*Kff*rs);
psrs = ps'*rs;
ws = (nrs^2*ps - (Kff*rs)*(psrs))/nrs^3;
NR = zeros(s, 1);
NR(s) = nrs;
PR = zeros(s, 1);
PR(s) = psrs;

qs = ws;
WV = zeros(s, s);
Kv = zeros(ndof, s);
for j = 1:(s - 1)
    vj = V(nf, j);
    wsvj = ws'*vj;
    Kvj = Kff*vj;
    qs = qs - Kvj*wsvj;
    WV(s, j) = wsvj;
    Kv(nf, j) = Kvj;
end
us = U(nf, s);
qsus = qs'*us;
nus = sqrt(us'*Kff*us);
zs = (R\(R'\(nus^2*qs - (Kff*us)*(qsus))))/nus^3;

NU = zeros(s, 1);
NU(s) = nus;
QU = zeros(s, 1);
QU(s) = qsus;

P(nf, s) = ps;
P(np, s) = zeros(size(np));
W(nf, s) = ws;
W(np, s) = zeros(size(np));
Q(nf, s) = qs;
Q(np, s) = zeros(size(np));
Z(nf, s) = zs;
Z(np, s) = zeros(size(np));
for i = (s-1):-1:1
    pi = -2*ff*(Vf(:, i)'*ff);
    vi = Vf(:, i);
    for j = (i + 1):s
        wj = W(nf, j);
        wjvi = wj'*vi;
        tj = T(nf, j);
        Ktj = Kff*tj;
        pi = pi - wjvi*Ktj - (vi'*Ktj)*wj;
        WV(j, i) = wjvi;
    end
    P(nf, i) = pi;
    
    ri = Rb(nf, i);
    nri = sqrt(ri'*Kff*ri);
    piri = pi'*ri;
    wi = (nri^2*pi - (Kff*ri)*(piri))/nri^3;
    W(nf, i) = wi;
    NR(i) = nri;
    PR(i) = piri;
    
    zip1 = Z(nf, i + 1);
    qi = wi - dKff*zip1;
    for j = 1:(i-1)
        vj = Vf(:, j);
        wivj = wi'*vj;
        qi = qi - Kv(nf, j)*wivj;
        WV(i, j) = wivj;
    end
    Q(nf, i) = qi;
    
    ui = U(nf, i);
    nui = sqrt(ui'*Kff*ui);
    qi = Q(nf, i);
    qiui = qi'*ui;
    NU(i) = nui;
    QU(i) = qiui;
    if i > 1
        zi = (R\(R'\(nui^2*qi - (Kff*ui)*qiui)))/nui^3;
        Z(nf, i) = zi;
    end
end

% Computing sensitivities
dc = zeros(1, nelm);
for elm = 1:nelm
    dofs = edof(elm, 2:end);
    ke = K0{elm};
    Zelm = Z(dofs, :);
    Uelm = U(dofs, :);
    Telm = T(dofs, :);
    Relm = Rb(dofs, :);
    Velm = V(dofs, :);
    
    dcelm = 0;
    for i = 1:s
        tie = Telm(:, i);
        if i < s
            dcelm = dcelm + Zelm(:, i + 1)'*ke*tie;
        end
        
        for j = 1:(i - 1)
            dcelm = dcelm + (tie'*ke*Velm(:, j))*WV(i, j); 
        end
        
        uie = Uelm(:, i);
        dcelm = dcelm + 1/2*QU(i)*(uie'*ke*uie)/NU(i)^3;
        
        rie = Relm(:, i);
        dcelm = dcelm + 1/2*PR(i)*(rie'*ke*rie)/NR(i)^3;
    end
    dc(elm) = dcelm;
end