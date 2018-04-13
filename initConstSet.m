function [ eps_r, ni, doping, phi, boundary, nn, pp, Dn, Dp ] ...
    = initConstSet( eps_r, ni, doping, phi, boundary, ...
                    nn, pp, Dn, Dp, Dn0, Dp0, ...
                    x_idx, x_int, x_dlt, z_idx, z_int, z_dlt, ...
                    Nd, Vt, Vss, Vds, Vgs )

ni0 = 1.0750038e+10;    % [cm^-3]
% ni0 = 1.0758720e+10;    % [cm^-3]
eps_r_si = 11.7;        % [1]
eps_r_ox = 3.9;         % [1]

%% index setting
% set index x-direction
nx = x_idx{end};
x_int1 = x_int{1};
x_int2 = x_int{2};

% set index z-direction
nz = z_idx{end};
z_int1 = z_int{1};
z_int2 = z_int{2};

% set metal line 
x_source = 1;   % left source 
z_source = [z_int1(1) z_idx{2} z_int2(2)];  

x_drain  = nx;  % right drain
z_drain  = z_source;              

x_gate1  = [x_int1(1) x_idx{2} x_int2(2)];  
z_gate1  = 1;   % lower gate 

x_gate2  = x_gate1;                         
z_gate2  = nz;  % upper gate

%% constant setting 
% constant setting : epsilon
eps_r(:,z_idx{1}) = eps_r_ox;    % [1] in oxide.
eps_r(:,z_idx{2}) = eps_r_si;    % [1] in semi.
eps_r(:,z_idx{3}) = eps_r_ox;    % [1] in oxide.

% constant setting : intrinsic carrier profile 
ni(:, z_idx{1}) = 0;             % [m^-3] in oxide 
ni(:, z_idx{2}(1:end)) = (ni0*1e+6);    % [m^-3] in semi.
% ni(:, z_idx{2}(1)) = 0;          % interface
% ni(:, z_idx{2}(end)) = 0;        % interface
ni(:, z_idx{3}) = 0;             % [m^-3] in oxide

% constant setting : doping profile 
doping(:, z_idx{1}) = 0;                % in oxide 
doping(x_idx{1},z_idx{2}(1:end)) = Nd(1)*1e+6;   % in semi. (n+)
doping(x_idx{2},z_idx{2}(1:end)) = Nd(2)*1e+6;   % in semi. (i)
doping(x_idx{3},z_idx{2}(1:end)) = Nd(3)*1e+6;   % in semi. (n+)
% doping(:, z_idx{2}(1)) = 0;             % interface 
% doping(:, z_idx{2}(end)) = 0;           % interface 
doping(:, z_idx{3}) = 0;                % in oxide
% 
doping(x_int2(1),z_idx{2}(1:end)) = Nd(1)*1e+6;   % in semi. (n+)
doping(x_int1(2),z_idx{2}(1:end)) = Nd(3)*1e+6;   % in semi. (n+)

% doping(:,z_idx{2}(1)) = 2*doping(:,z_idx{2}(1));   % in semi. (n+)
% doping(:,z_idx{2}(end)) = 2*doping(:,z_idx{2}(end));   % in semi. (n+)

% doping(:,z_idx{2}(1)) = 0;
% doping(:,z_idx{2}(end)) = 0;

% Dn profile  
helper = zeros(nx,nz);
helper(:, z_idx{1}) = 0;        % in oxide 
helper(:, z_idx{2}) = Dn0;      % in semi.
helper(:, z_idx{3}) = 0;        % in oxide 

Dn(:, :, 1) = helper;     % left
Dn(:, z_idx{2}(1), 1)   = Dn0*z_dlt(z_int2(1))/( z_dlt(z_int1(1)-1) + z_dlt(z_int2(1)) );
Dn(:, z_idx{2}(end), 1) = Dn0*z_dlt(z_int1(2)-1)/( z_dlt(z_int1(2)-1) + z_dlt(z_int2(2)) );

Dn(:, :, 2) = helper;     % right
Dn(:, z_idx{2}(1), 2)   = Dn0*z_dlt(z_int2(1))/( z_dlt(z_int1(1)-1) + z_dlt(z_int2(1)) );
Dn(:, z_idx{2}(end), 2) = Dn0*z_dlt(z_int1(2)-1)/( z_dlt(z_int1(2)-1) + z_dlt(z_int2(2)) );

Dn(:, :, 3) = helper;     % lower
Dn(:, z_idx{2}(1), 3)   = 0;

Dn(:, :, 4) = helper;     % upper
Dn(:, z_idx{2}(end), 4) = 0;

% Dp profile 
helper = zeros(nx,nz);
helper(:, z_idx{1}) = 0;        % in oxide 
helper(:, z_idx{2}) = Dp0;      % in semi.
helper(:, z_idx{3}) = 0;        % in oxide 

Dp(:, :, 1) = helper;     % left
Dp(:, z_idx{2}(1), 1)   = Dp0*z_dlt(z_int2(1))/( z_dlt(z_int1(1)-1) + z_dlt(z_int2(1)) );
Dp(:, z_idx{2}(end), 1) = Dp0*z_dlt(z_int1(2)-1)/( z_dlt(z_int1(2)-1) + z_dlt(z_int2(2)) );

Dp(:, :, 2) = helper;     % right
Dp(:, z_idx{2}(1), 2)   = Dp0*z_dlt(z_int2(1))/( z_dlt(z_int1(1)-1) + z_dlt(z_int2(1)) );
Dp(:, z_idx{2}(end), 2) = Dp0*z_dlt(z_int1(2)-1)/( z_dlt(z_int1(2)-1) + z_dlt(z_int2(2)) );

Dp(:, :, 3) = helper;     % lower
Dp(:, z_idx{2}(1), 3)   = 0;

Dp(:, :, 4) = helper;     % upper
Dp(:, z_idx{2}(end), 4) = 0;

% Dirichlet boundary condition (phi)
boundary(x_gate1, z_gate1, 1) = Vgs;    % gate (lower) 
boundary(x_gate2, z_gate2, 1) = Vgs;    % gate (upper)
% boundary(x_source,z_source,1) = Vss;    % source (left)
% boundary(x_drain, z_drain, 1) = Vss;    % drain (right)

% initial phi guess 
phi = initPotential(phi, doping, ni, Vt);

BC_index = find(boundary ~= 0); 
phi(BC_index) = boundary(BC_index);

% initial electron / hole setting 
% p-type
pp_index = find(doping < 0); 
pp(pp_index) = -doping(pp_index) + ni(pp_index);
nn(pp_index) = (ni0*1e+6).^2./pp(pp_index);
% n-type
nn_index = find(doping > 0);
nn(nn_index) = doping(nn_index) + ni(nn_index);
pp(nn_index) = (ni0*1e+6).^2./nn(nn_index);
% i-type
ii_index = find(doping == 0);
nn(ii_index) = ni(ii_index);
pp(ii_index) = ni(ii_index);

% % boundary condiiton (source)
% nn_s = (sqrt((Nd(1)*1e+6)^2 + 4*(ni0*1e+6)^2) + (Nd(1)*1e+6)) / 2;
% pp_s = (ni0*1e+6)^2/nn_s;
% boundary(x_source, z_source, 2) = nn_s;
% boundary(x_source, z_source, 3) = pp_s;
% 
% % boundary condiiton (drain)
% nn_d = (sqrt((Nd(3)*1e+6)^2 + 4*(ni0*1e+6)^2) + (Nd(3)*1e+6)) / 2;
% pp_d = (ni0*1e+6)^2/nn_d;
% boundary(x_drain, z_drain, 2) = nn_d;
% boundary(x_drain, z_drain, 3) = pp_d;

end

