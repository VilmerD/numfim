function y = gsorth(x, V, A)
%GSORTH computes the gram-schmidt orthogonalization of x with respect to
% the vectors in V in the A-energy norm

if nargin < 3
    A = 1;
end

y = x;
for k = 1:size(V, 2)
    vk = V(:, k);
    ak = (x'*A*vk);
    y = y - ak*vk;
end

end