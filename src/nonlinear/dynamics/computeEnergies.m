function [Kin, Int, Intc] = computeEnergies(sys, rho, u, du)
% Computes the energies in the system
%
% Inputs:
%   sys:    init2D object
%   rho:    density of system
%   u:      displacement field  [m x 1]
%   du:     velocity field      [m x 1]

[~, n] = size(u);
Kin = zeros(1, n); Int = Kin; Intc = Kin;

wbar = waitbar(0, '1', 'Name', 'Computing energies', 'CreateCancelBtn', ...
        'setappdata(gcbf, ''canceling'', 1)');
setappdata(wbar, 'canceling', 0);

for k = 1:n
    if getappdata(wbar, 'canceling')
        break
    end
   [Ke, Ie] = sys.energy(rho, u(:, k), du(:, k));
   Kin(k) = Ke; Int(k) = Ie; Intc(k) = sys.energyc(u(:, k));
   waitbar(k/n, wbar, 'Please wait.')
end
delete(wbar)
end