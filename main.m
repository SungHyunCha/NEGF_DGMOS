clear;

%% sequence 
% #0. Define basic information (mesh size / position / constant) 
% #1. Solve 2D init Poisson eq.
% #2. Repeat following loop
%   #2-1. Solve Poisson eq, Current eq (electron, hole)
%   #2-2. If error is larger than 1e-5, repeat #2. 
% #3. Print result such as concentration, potential

do_init = 0;

if do_init == 1
    %% model parameters 
    Vg_bias = 0.0; 
%     Vgate  = 0.533744;
    Vgate  = 0.5100;
    Vsource = 0.6113137; 
    Vgs = Vgate + Vg_bias; 
    Vss = Vsource;
    Vds = Vsource;
    
    Nd = [2e+20 ; 0 ; 2e+20]; 

    %% position grid setting 
    x = csvread('xmesh_0.5.csv', 1, 0, [1 0 57 0])*1e+3;   % node position [nm]
    x_nod = size(x,1);
    x_dlt = (x(2:end) - x(1:end-1));  % node spacing [nm]
    x_int1 = find(x_dlt == 0);        % at least two interface 
    x_int2 = x_int1 + 1;
    x_int = {x_int1 x_int2};

    x_idx = {1:x_int1(1)};            % cell notation 
    for i = 2:size(x_int1,1)
        x_idx{end+1} = x_int2(i-1):x_int1(i);
    end
    x_idx{end+1} = x_int2(i):x_nod;
    x_idx{end+1} = x_nod;

    z = csvread('zmesh_0.125.csv', 1, 0, [1 0 59 0])*1e+3;
    z_nod = size(z,1);
    z_dlt = (z(2:end) - z(1:end-1));
%     z_dlt = round(z_dlt*1e+3)*1e-3;
    z_int1 = find(z_dlt == 0);
    z_int2 = z_int1 + 1;
    z_int = {z_int1 z_int2};

    z_idx = {1:z_int1(1)};
    for i = 2:size(z_int1,1)
        z_idx{end+1} = z_int2(i-1):z_int1(i);
    end
    z_idx{end+1} = z_int2(i):z_nod;
    z_idx{end+1} = z_nod;

    % delta setting at interface 
    x_dlt(x_int1) = x_dlt(x_int2);
    z_dlt(z_int1) = z_dlt(z_int2);

    %% constant (position independent)
    q = 1.602192e-19;       % [J/eV] or [C]
    Vt = 1.380662e-23*300/q;% [V]
    eps0 = 8.8542e-12;  	% [F/m]
    Egap = 1.11;            % [eV] band gap

    %% constant setting 
    % Constant Setting (Poisson eq.)
    % related with Poisson eq. 
    eps_r      = zeros(x_nod,z_nod);   % relative permitivity (left, right, bottom, top)
    ni         = zeros(x_nod,z_nod);     % intrinsic electron concentration 
    doping     = zeros(x_nod,z_nod);     % doping concentration [m^-3]
    boundary   = zeros(x_nod,z_nod,3);   % Dirichlet boundary (phi, n, p)
    phi        = zeros(x_nod,z_nod);     % potential [V] 

    % Constant Setting (2D continuity eq.)
    nn = zeros(x_nod,z_nod);         % electron concentration [m^-3]
    pp = zeros(x_nod,z_nod);         % hole     concentration [m^-3]
    Dn = zeros(x_nod,z_nod,4);
    Dp = zeros(x_nod,z_nod,4);
    Dn0 = 1417e-4*Vt;       % electron diff. coefficient [m^2/sec]
    Dp0 = 470.5e-4*Vt;      % hole     diff. coefficient [m^2/sec]

    % Constant Setting 
    [ eps_r, ni, doping, phi, boundary, nn, pp, Dn, Dp ] ...
        = initConstSet( eps_r, ni, doping, phi, boundary, ...
                        nn, pp, Dn, Dp, Dn0, Dp0, ...
                        x_idx, x_int, x_dlt, z_idx, z_int, z_dlt, ...
                        Nd, Vt, Vss, Vds, Vgs );
                    
    %% initial Poisson eq. 
    jbase = configueJbase(eps_r, boundary, x_idx, x_int, x_dlt, z_idx, z_int, z_dlt);

    [phi, nn, pp ] = initPoisson2D( 100, jbase, ni, phi, boundary(:,:,1), doping, ...
                          x_dlt, x_int, z_dlt, z_int, eps0, Vt, q);
    save('init.mat');
    return;
else 
    load('init.mat');
end
% bias = 0.6;
% nVds = 7;   
% Vds = linspace(0, bias, nVds);     % simulation target 
bias = 0.6;
nVgs = 4;   
Vgs = linspace(0, bias, nVgs);     % simulation target 

nn_new = nn;
phi_new = phi;

for i = 1:1
    for j = 1:100
        phi_old = phi_new;
        % delta1, nodenum2, valley
        tic;
%         n1 = 2*schOneValley(0.2e-3, 1001, 1, phi_new(:,z_idx{2}), Vds(i), z(z_idx{2}), x, x_int2, Egap, Vt, q); 
%         n2 = 2*schOneValley(0.2e-3, 1001, 2, phi_new(:,z_idx{2}), Vds(i), z(z_idx{2}), x, x_int2, Egap, Vt, q); 
%         n3 = 2*schOneValley(0.2e-3, 1001, 3, phi_new(:,z_idx{2}), Vds(i), z(z_idx{2}), x, x_int2, Egap, Vt, q); 
        n1 = 2*schOneValley(0.2e-3, 1001, 1, phi_new(:,z_idx{2}), 0.0, z(z_idx{2}), x, x_int2, Egap, Vt, q); 
        n2 = 2*schOneValley(0.2e-3, 1001, 2, phi_new(:,z_idx{2}), 0.0, z(z_idx{2}), x, x_int2, Egap, Vt, q); 
        n3 = 2*schOneValley(0.2e-3, 1001, 3, phi_new(:,z_idx{2}), 0.0, z(z_idx{2}), x, x_int2, Egap, Vt, q); 


        nn_new(:,z_idx{2}) = n1 + n2 + n3;
%         save('result.mat');
%         elapsedTime = toc
%         return; 

%         [ phi_new ] = nLinPoisson2D( 100, jbase, nn_new, ni, phi_new, boundary(:,:,1), doping, ...
%                               x_dlt, x_int, z_dlt, z_int, eps0, Vt, q);
        [ phi_new ] = nLinPoisson2D( 100, jbase, nn_new, ni, phi_new, boundary(:,:,1), Vgs(i), doping, ...
                              x_dlt, x_int, z_dlt, z_int, eps0, Vt, q);

        elapsedTime(i,j) = toc;
%         disp(sprintf('Elapsed Time(sec): %03d \n', elapsedTime(i,j)));
        
        stop(i,j) = max(max(abs(phi_new-phi_old)));
        disp(sprintf('[%d]self-consist loop[%d]-error: %d \n', i, j, stop(i,j)));
        if (stop(i,j) < 1e-4) || (j == 100)
            save(sprintf('result_%03d.mat',i));
            break;
        end
    end 
end 
    
save('result.mat');

return;

nn = nn_new;
phi = phi_new;

ii = 26;
% potential comparison plot 
figure('Color', [1,1,1]); % z-slice
plot(x, phi_new(:,ii), 'b'); hold on; plot(x, phi(:,ii), 'r');
xlabel('Position (nm)'); ylabel('Potential (V)'); 
title('Electric potential slice(x-direction) profile'); grid on;
legend_h = legend('1mV 10step', '10mV 1step');
set(legend_h,'Location','Best','FontSize',11);

% error plot! 
figure('Color', [1,1,1]); % z-slice
plot(x, 100*abs( (phi_new(:,ii) - phi(:,ii) )./phi_new(:,ii) ), '--ro');
xlabel('Position (nm)'); ylabel('Potential ratio (%)'); 
title('Potential ratio profile'); grid on;
legend_h = legend('abs((A-B)/A)');
set(legend_h,'Location','Best','FontSize',11);
% find maximum error! 
max(100*abs( (phi_new(:,ii) - phi(:,ii) )./phi_new(:,ii)))

ii = 26;
% concentration comparison plot 
figure('Color', [1,1,1]); % z-slice
% plot(x, nn_new(:,ii)*1e-6, 'b'); hold on; plot(x, nn(:,ii)*1e-6, 'r');
plot(x, nn_new(:,ii)/max(nn_new(:,ii)), 'b'); hold on; plot(x, nn(:,ii)/max(nn(:,ii)), 'r');
xlabel('Position (nm)'); ylabel('Electron concentration (cm^-3)'); 
title('Electron concentration profile'); grid on;
legend_h = legend('1mV 10step', '10mV 1step');
set(legend_h,'Location','Best','FontSize',11);

% error plot! 
figure('Color', [1,1,1]); % z-slice
plot(x, 100*abs( (nn_new(:,ii) - nn(:,ii) )./nn_new(:,ii) ), '--ro');
xlabel('Position (nm)'); ylabel('Electron concentration ratio (%)'); 
title('Electron concentration ratio profile'); grid on;
legend_h = legend('abs((A-B)/A)');
set(legend_h,'Location','Best','FontSize',11);
% find maximum error! 
max(100*abs( (nn_new(:,ii) - nn(:,ii) )./nn_new(:,ii) ))




