function [ integration ] = Fhalfinv( E, nE, hBar, delta, m, Vt, q)
tOption = 0; % transpose option 
if size(E,2) == 1   % column vector 
    tOption = 0;
elseif size(E,1) == 1   % row vector 
    tOption = 1;
    E = E';
else
    disp('Error: Invalid Input Format (Fermi integral)');
    integration = 0;
    return;
end

eps = linspace(0, 10, nE);
d_eps = eps(2) - eps(1);
eps = eps(1:end-1) + d_eps/2;
[eps_mat, E_mat] = meshgrid(eps, E);
eps2_mat = eps_mat + d_eps/2;
eps1_mat = eps_mat - d_eps/2;

integration = (hBar*delta)^-1*sqrt(m*q*Vt/2)/pi*sum( 2*(1./(1+exp(eps_mat - E_mat/Vt))).*(sqrt(eps2_mat) - sqrt(eps1_mat)), 2 );

if tOption == 1 % recover transposed array 
    integration = integration';
end 
end
