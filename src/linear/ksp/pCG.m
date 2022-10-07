function [x, nits] = pCG(Afun, b, C, rtol, x0)
% CG finds the vector xn which minimizes the error en in the A-norm
%

if nargin < 4
    rtol = 1e-9;
end

if nargin < 5
    x0 = zeros(size(b));
end

x = C*x0;
afunp = @(x) C'\(Afun((C\x)));
r = C'\(b - Afun(x0));
p = r;

rm1trm1 = r'*r;
nits = 0;
while rm1trm1 > rtol
    Ap = afunp(p);
    alph = rm1trm1/(p'*Ap);
    x = x + alph*p;
    r = r - alph*Ap;
    rtr = r'*r;
    bet = rtr/rm1trm1;
    p = r + bet*p;
    
    rm1trm1 = rtr;
    nits = nits + 1;
end
x = C\x;
end