function [x, nits] = CG(Afun, b, x0)
% CG finds the vector xn which minimizes the error en in the A-norm
%

if nargin == 2
    x0 = zeros(size(b));
end

x = x0;
r = b - Afun(x);
p = r;

rm1trm1 = r'*r;
nits = 0;
while rm1trm1 > 1e-6
    Ap = Afun(p);
    alph = rm1trm1/(p'*Ap);
    x = x + alph*p;
    r = r - alph*Ap;
    rtr = r'*r;
    bet = rtr/rm1trm1;
    p = r + bet*p;
    
    rm1trm1 = rtr;
    nits = nits + 1;
end
end