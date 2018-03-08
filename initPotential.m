function [ phi ] = initPotential( phi, doping, ni, Vt )

pp_index = find(doping < 0);    % p-type
phi(pp_index) = -(Vt)*log(-doping(pp_index)./ni(pp_index));

nn_index = find(doping > 0);    % n-type
phi(nn_index) = +(Vt)*log(+doping(nn_index)./ni(nn_index));
end

