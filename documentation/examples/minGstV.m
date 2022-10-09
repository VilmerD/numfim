%% Minimum max gamma (1/lambda) s.t. volume
% Approximation of max gamma using p-norm 
% Single density field with eta=0.5 

% OPEN ISSUES
% EFFICIENCY OF LOOPS 
% BCS FOR LOADS AND SUPPORTS
% INTERPOLATION SCHEMES?
% ADD CA BASED IN VILMER'S CODE


clear;
close all;
commandwindow;
clc;
tic;
% prtype = 'WEIGHTED';
% weight = 1;
prtype = 'BLF'; 
% prtype =  'COMP';

%% Domain size and discretization
domain = 'column';      % other options are: 'ubracket' 'biclamped' 
sizex = 200;            % physical size in x-direction
sizey = 40;             % physical size in y-direction
helem = 1;              % element size (all elements are square)
%% Optimization parameters
pE = 3;                 % SIMP penalty for linear stiffness
pS = 3;                 % SIMP penalty for stress stiffness
pN = 8;                 % p-norm / KS exponent
rmin = 3.1;             % filter radius (for convolution density filter)
nevals = 6;				% number of eigenvalues to consider
filename = 'column_trial.mat';
%% Material properties
Emax = 2e5;
Emin = Emax*1e-6;
nu = 0.3;
%% Prepare design domain
% X = nodal coordinates
% T = element connectivity and data
% i_img,j_img = for displaying design as a matrix using imagesc
switch domain
    case 'column'
        [X,T,i_img,j_img,solids,voids,F,freedofs] = generate_column(sizex,sizey,helem,false);
    case 'biclamped'
        error('bi-clamped domain not supported yet');
        % [X,T,i_img,j_img] = generate_biclamped(sizex,sizey,helem,false);
    case 'lbracket'
        % [X,T,i_img,j_img] = generate_lbracket(sizex,sizey,helem,false);
         error('l-bracket domain not supported yet');
end
%% Prepare FEA (88-line style)
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nelem = size(T,1);
nodes = size(X,1);
ndof = 2*nodes;
edofMat = zeros(nelem,8);
edofMat(:,1:2) = [2*T(:,1)-1 2*T(:,1)];
edofMat(:,3:4) = [2*T(:,2)-1 2*T(:,2)];
edofMat(:,5:6) = [2*T(:,3)-1 2*T(:,3)];
edofMat(:,7:8) = [2*T(:,4)-1 2*T(:,4)];
iK = reshape(kron(edofMat,ones(8,1))',64*nelem,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelem,1);
iF = reshape(edofMat',8*nelem,1);
U = zeros(ndof,1); PHI = zeros(ndof,nevals);
% For stress computations 
B = 1/helem*[-1/2 0 1/2 0 1/2 0 -1/2 0
    0 -1/2 0 -1/2 0 1/2 0 1/2
    -1/2 -1/2 -1/2 1/2 1/2 1/2 1/2 -1/2];           % strain-displacement matrix
D = 1/(1-nu^2)*[ 1 nu 0;nu 1 0;0 0 (1-nu)/2];       % constitutive matrix - plane stress
% For stress stifffness matrix
BNL = 1/helem*[-1/2 0 1/2 0 1/2 0 -1/2 0
               -1/2 0 -1/2 0 1/2 0 1/2 0
               0 -1/2 0 1/2 0 1/2 0 -1/2
               0 -1/2 0 -1/2 0 1/2 0 1/2];          % BNL matrix
%% Prepare augmented PDE filter (with Robin BC) 
xi = 1;                     % default for Robin BC
xi_corner = 1;             	% large xi near re-entrant corner
l_o = rmin/2/sqrt(3);       % bulk length scale parameter
l_s = rmin/2/sqrt(3)*xi;    % default surface length scale parameter
% PDE filter stiffness matrix 
edofMatF = [(edofMat(:,7)+1)/2 (edofMat(:,5)+1)/2 (edofMat(:,3)+1)/2 (edofMat(:,1)+1)/2];
KEF = l_o^2*[4 -1 -2 -1; -1  4 -1 -2; -2 -1  4 -1; -1 -2 -1  4]/6 + ...
 [4  2  1  2;  2  4  2  1;  1  2  4  2;  2  1  2  4]/36; 
iKF = reshape(kron(edofMatF,ones(4,1))',16*nelem,1);
jKF = reshape(kron(edofMatF,ones(1,4))',16*nelem,1);
sKF = reshape(KEF(:)*ones(1,nelem),16*nelem,1);
% PDE filter "mass" matrix
sMF = 0*sKF;
ME1 = l_s*[2 0 0 1; 0 0 0 0; 0 0 0 0; 1 0 0 2]/6;                        % mass matrix left face
ME2 = l_s*[0 0 0 0; 0 0 0 0; 0 0 2 1; 0 0 1 2]/6;                        % mass matrix bottom face
ME31 = l_s*[0 0 0 0; 0 2 1 0; 0 1 2 0; 0 0 0 0]/6;                       % mass matrix right face, standard
ME32 = l_s*(xi_corner/xi)*[0 0 0 0; 0 2 1 0; 0 1 2 0; 0 0 0 0]/6;        % mass matrix right face, corner
ME41 = l_s*[2 1 0 0; 1 2 0 0; 0 0 0 0; 0 0 0 0]/6;                       % mass matrix top face, standard
ME42 = l_s*(xi_corner/xi)*[2 1 0 0; 1 2 0 0; 0 0 0 0; 0 0 0 0]/6;        % mass matrix top face, corner
% Populate "mass" matrix  
idv = zeros(1,nelem);                               % identifies elements for the boundary surface terms
idv = idv*0; idv(T(:,9)==1) = 1;                    % left face
sMF = sMF + reshape(ME1(:)*idv,16*nelem,1);
idv = idv*0; idv(T(:,10)==1) = 1;                   % bottom face
sMF = sMF + reshape(ME2(:)*idv,16*nelem,1);
idv = idv*0; idv(T(:,11)==1) = 1; idv(solids)=0;    % right face standard
sMF = sMF + reshape(ME31(:)*idv,16*nelem,1);
idv = idv*0; idv(T(:,11)==2) = 1; idv(solids)=0;    % right face corner 
sMF = sMF + reshape(ME32(:)*idv,16*nelem,1);
idv = idv*0; idv(T(:,12)==1) = 1; idv(solids)=0;    % top face standard 
sMF = sMF + reshape(ME41(:)*idv,16*nelem,1);
idv = idv*0; idv(T(:,12)==2) = 1; idv(solids)=0;    % top face corner 
sMF = sMF + reshape(ME42(:)*idv,16*nelem,1);
KF = sparse(iKF,jKF,sKF+sMF);
LF = chol(KF,'lower');
% Transformation Matrix
iTF = reshape(edofMatF,4*nelem,1);
jTF = reshape(repmat(1:nelem,4,1)',4*nelem,1);
sTF = repmat(1/4,4*nelem,1);
TF = sparse(iTF,jTF,sTF);
%% Initialize optimization
maxloop = 200;
volfrac = 0.5;
x = volfrac*ones(nelem,1);
x(T(:,5)==1) = 1e-6;    % voids
x(T(:,5)==2) = 1-1e-6;  % solids
beta = 1;
njumps = 8;
betamax = 16;
dbeta = (betamax/beta)^(1/njumps); % factor for multiplying beta
pace = min(20,maxloop/(njumps+1));
loop = 0;
Stats = zeros(maxloop,20);
%% Initialize MMA
m     = 1;                  % number of general constraints.
n     = nelem;              % number of design variables x_j.
xmin  = 1e-6*ones(n,1);     % column vector with the lower bounds for the variables x_j.
xmin(solids) = 1-1e-3;      % lower bound for solids
xmax  = ones(n,1);          % olumn vector with the upper bounds for the variables x_j.
xmax(voids) = 1e-3;         % upper bound for voids
xold1 = x(:);               % xval, one iteration ago (provided that iter>1).
xold2 = x(:);               % xval, two iterations ago (provided that iter>2).
low   = 0*ones(n,1);        % column vector with the lower asymptotes from the previous iteration (provided that iter>1).
upp   = ones(n,1);          % column vector with the upper asymptotes from the previous iteration (provided that iter>1).
a0    = 1;                  % the constants a_0 in the term a_0*z.
a     = zeros(m,1);         % column vector with the constants a_i in the terms a_i*z.
c_MMA = 1000*ones(m,1);     % column vector with the constants c_i in the terms c_i*y_i.
d     = ones(m,1);          % column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
fval = ones(m,1);           % initial constraint values
%% Start iteration 
while ((loop < maxloop || beta<0.99*betamax) || fval(1,1) > 0)
    loop = loop + 1;
    %% PDE filtering and single projection
    xTilde = (TF'*(LF'\(LF\(TF*x(:)))));
    xPhys = (tanh(beta*0.5)+tanh(beta*(xTilde-0.5)))/...
        (tanh(beta*0.5)+tanh(beta*(1-0.5)));   
    %% Solve static equation
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^pE*(Emax-Emin)),64*nelem,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    R = chol(K(freedofs,freedofs));
    U(freedofs) = R\(R'\F(freedofs));
    % SOLVE STATICS WITH CA
    %
    %
    %
    
    %% Compliance and its sensitivity
    ce = sum((U(edofMat)*KE).*U(edofMat),2);
    comp = sum(sum((Emin+xPhys.^pE*(Emax-Emin)).*ce));
    dc = -pE*(Emax-Emin)*xPhys.^(pE-1).*ce;
    dv = ones(nelem,1);
    %% Compute geometric stiffness 
    EPS = U(edofMat)*B';                                % strain
    SIG = EPS*D;                                        % stress for E=1
    sG = 0*sK; dsG = sG;
    for el = 1:nelem
        tau = [SIG(el,1) SIG(el,3) 0 0;
            SIG(el,3) SIG(el,2) 0 0 ;
            0 0 SIG(el,1) SIG(el,3);
            0 0 SIG(el,3) SIG(el,2)];
        BNLTtau = BNL'*tau;
        BNLTtauBNL = BNLTtau*BNL;
        l1 = (el-1)*64+1; l2 = el*64;
        sG(l1:l2) = Emax*xPhys(el,1)^pS*helem^2*BNLTtauBNL(:);
		dsGdx(l1:l2) = Emax*pS*xPhys(el,1)^(pS-1)*helem^2*BNLTtauBNL(:);
    end
    %% Solve eigenvalue problem
    KNL = sparse(iK,jK,sG); KNL = (KNL+KNL')/2;
    Kfree = K(freedofs,freedofs); 
    KNLfree = KNL(freedofs,freedofs); 
    [evecs,evals] = eigs(KNLfree,-Kfree,nevals,'la');
    % SOLVE EIGS WITH CA
    %
    %
    %
    
    mu = diag(evals);
    mu_max_acc = max(mu);
    mu_max_app = norm(mu,pN); 
    lambda = 1./mu;
    lambda_min_acc = 1/mu_max_acc;
    lambda_min_app = 1/mu_max_app;
    % Enforce orthonormality
    phiKphi = evecs'*(Kfree*evecs);
	PHI(freedofs,:) = evecs./sqrt(diag(phiKphi)');
	%% Sensitivities of buckling load factor
    % Compute sensitivity term #1 -PHI'*(dGdx+mu*dKdx)*PHI
	dmu1 = zeros(nelem,nevals);
    for el = 1:nelem
        l1 = (el-1)*64+1; l2 = el*64;
        dsGe = reshape(dsGdx(l1:l2),8,8);
        dsKe = pE*(Emax-Emin)*xPhys(el,1)^(pE-1)*KE;
        dmu11 = sum(PHI(edofMat(el,:)',:).*(dsGe*PHI(edofMat(el,:)',:)));
        dmu12 = mu'.*sum(PHI(edofMat(el,:)',:).*(dsKe*PHI(edofMat(el,:)',:)));
        dmu1(el,:) = - (dmu11 + dmu12); 
    end
	% Compute adjoint loads
	dSIGdu = D*B; % 3*8
    dSIGxxdu = dSIGdu(1,:);
    dSIGyydu = dSIGdu(2,:);
    dSIGxydu = dSIGdu(3,:);
    ADJload = 0*PHI;
    for el = 1:nelem
        edof = edofMat(el,:);
        for i = 1:8 % prepare 8*8 dsGdu separately for each DOF
            dtaudu = [dSIGxxdu(i) dSIGxydu(i) 0 0;
                dSIGxydu(i) dSIGyydu(i) 0 0 ;
                0 0 dSIGxxdu(i) dSIGxydu(i);
                0 0 dSIGxydu(i) dSIGyydu(i)];
            dsGdu = Emax*xPhys(el,1)^pS*helem^2*(BNL'*(dtaudu*BNL)); % 8*8  
            dsGduPHI = dsGdu*PHI(edof,:); % 8*nevals
            PHIdsGduPHI = sum(PHI(edof,:).*dsGduPHI); % 1*nevals
            % This gives the contribution of element el to a specific u
            ind = (el-1)*8+i;
            ADJload(edof(i),:) = ADJload(edof(i),:) - PHIdsGduPHI;
        end
    end
    ADJsol = 0*ADJload;
    ADJsol(freedofs,1:nevals) =  R\(R'\ADJload(freedofs,1:nevals));
    % SOLVE ADJOINT WITH CA
    %
    %
    %
    
    dmu2 = 0*dmu1;
    for j = 1:nevals
        adjsol = ADJsol(:,j);
        vals = sum((adjsol(edofMat)*KE).*U(edofMat),2);
        dmu2(:,j) = -pE*(Emax-Emin)*xPhys.^(pE-1).*vals;
    end
	dmu = dmu1 + dmu2;
    term1 = 1/pN*(sum(mu.^pN))^(1/pN-1);
    term2 = (mu.^(pN-1))'*dmu';
    dmu_max_app = term1*pN*term2';
     %% Chain rule for projection and filter
    dxphys = (1 - (tanh(beta*(xTilde(:)-0.5))).^2)*beta / ...
        (tanh(beta*0.5)+tanh(beta*(1-0.5)));
    dv(:) = dv(:).*dxphys(:);
    dc(:) = dc(:).*dxphys(:);
    dmu_max_app = dmu_max_app(:).*dxphys(:);
    dv(:) = TF'*(LF'\(LF\(TF*dv(:))));
    dc(:) = TF'*(LF'\(LF\(TF*dc(:))));
    dmu_max_app(:) = TF'*(LF'\(LF\(TF*dmu_max_app(:))));
    %% Draw design and stress
    figure(1);
    clf;
    v_img = xPhys;
    top_img = sparse(i_img,j_img,v_img);
    subplot(3,1,1);
    imagesc(top_img);
    axis equal;
    axis tight;
    axis off;
    title('xPhys');
    v_img = dmu_max_app;
    top_img = sparse(i_img,j_img,v_img);
    subplot(3,1,2);
    imagesc(top_img);
    axis equal;
    axis tight;
    axis off;
    title('d{\mu}_{max}');
    v_img = dc;
    top_img = sparse(i_img,j_img,v_img);
    subplot(3,1,3);
    imagesc(top_img);
    axis equal;
    axis tight;
    axis off;
    title('d{c}');
    drawnow;
    %% MMA
    xval  = x;
    switch prtype
        case 'COMP' % minCstV
            f0val = comp;
            if (loop==1)
                scale = 10/f0val;
            end
            df0dx = scale*dc(:);
        case 'BLF' % min gamma
            f0val = mu_max_app;
            if (loop==1)
                scale = 10/f0val;
            end
            df0dx = scale*dmu_max_app(:);
        case 'WEIGHTED'
            f0val = weight*mu_max_app + (1-weight)*comp;
            if (loop==1)
                scale1 = 10/mu_max_app;
                scale2 = 10/comp;
            end
             df0dx = scale1*weight*dmu_max_app(:) + scale2*(1-weight)*dc(:);
    end
    fval(1,1) = mean(xPhys)/volfrac - 1;
    dfdx(1,:) = dv(:)/volfrac/n;
    [xmma,~,~,~,~,~,~,~,~,low,upp] = ...
        mmasub(m,n,loop,xval,max(xmin,xval-0.2),min(xmax,xval+0.2),xold1,xold2, ...
        f0val,df0dx,fval,dfdx,low,upp,a0,a,c_MMA,d); 
    % Update MMA Variables
    xnew     = xmma;
    xold2    = xold1(:);
    xold1    = xval(:);
    change = max(abs(xnew(:)-xval(:)));
    x = xnew;
    %% Continuation
    if (mod(loop,pace) == 0)
        beta = min(beta*dbeta,betamax);
    end
    %% Print results
    fprintf(' ITER: %3i BLF_acc: %6.3e BLF_app: %6.3e OBJ: %6.3e CONST: %6.3e CH: %4.3f BETAHS: %6.3f\n',...
        loop,lambda_min_acc,lambda_min_app,f0val,fval(1,1),change,beta);
    %% Save data
    Stats(loop,1:11) = [lambda_min_acc lambda_min_app f0val fval' beta lambda'];
end
runtime = toc; 
save(filename);



