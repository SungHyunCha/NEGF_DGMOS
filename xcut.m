% %% error check (without node interpolation)
% z = csvread('zmesh.csv', 1, 0, [1 0 21 0])*1e+3;
% z_nod = size(z,1);
% 
% % csv value rearrange! 
% ii = 35;
% % csv_phi = csvread('Vds_phi_0.0V_xcut.csv',1,0,[1 0 z_nod 1]);
% csv_phi = csvread('phi_xcut_Vgs_0.6V_Vds_0.05V.csv',1,0,[1 0 z_nod 1]);
% csv_nn  = csvread('nn_xcut_Vgs_0.6V_Vds_0.05V.csv',1,0,[1 0 z_nod 1]);
% csv_pp  = csvread('pp_xcut_Vgs_0.6V_Vds_0.05V.csv',1,0,[1 0 z_nod 1]);
% csv_mob = csvread('emob_xcut_Vgs_0.6V_Vds_0.05V.csv',1,0,[1 0 z_nod 1]);

%% phi plot
ii = 10;    % 10, 29, 48

nn = nn_new;
phi = phi_new;

% notice: '??' should be changed! 
figure('Color', [1,1,1]); % z-slice
plot(z, phi(ii,:), 'b');
xlabel('Position (nm)'); ylabel('Potential (V)'); 
title('Electric potential slice(z-direction) profile'); grid on;
% 
% return;

figure('Color', [1,1,1]); % z-slice
plot(z, nn(ii,:)'*1e-6, 'b');
xlabel('Position (nm)'); ylabel('Electron concentration (cm^-3)'); 
title('Electron concentration profile'); grid on;
return; 

%% nn plot
% notice: '??' should be changed! 
figure('Color', [1,1,1]); % z-slice
plot(z(6:46), n1(ii,:)'*1e-6, 'b');
xlabel('Position (nm)'); ylabel('Electron concentration (cm^-3)'); 
title('Electron concentration profile'); grid on;

figure('Color', [1,1,1]); % z-slice
plot(z(6:46), nn1(ii,:)'*1e-6, 'b');
xlabel('Position (nm)'); ylabel('Electron concentration (cm^-3)'); 
title('Electron concentration profile'); grid on;

figure('Color', [1,1,1]); % z-slice
plot(z(6:46), n2(ii,:)'*1e-6, 'b');
xlabel('Position (nm)'); ylabel('Electron concentration (cm^-3)'); 
title('Electron concentration profile'); grid on;

figure('Color', [1,1,1]); % z-slice
plot(z(6:46), nn2(ii,:)'*1e-6, 'b');
xlabel('Position (nm)'); ylabel('Electron concentration (cm^-3)'); 
title('Electron concentration profile'); grid on;

figure('Color', [1,1,1]); % z-slice
plot(z(6:46), n3(ii,:)'*1e-6, 'b');
xlabel('Position (nm)'); ylabel('Electron concentration (cm^-3)'); 
title('Electron concentration profile'); grid on;

figure('Color', [1,1,1]); % z-slice
plot(z(6:46), nn3(ii,:)'*1e-6, 'b');
xlabel('Position (nm)'); ylabel('Electron concentration (cm^-3)'); 
title('Electron concentration profile'); grid on;


return;

% %% pp plot
% qx = csv_pp(:,1)*1e+3;
% qpp = csv_pp(:,2);
% % notice: '??' should be changed! 
% figure('Color', [1,1,1]); % z-slice
% plot(qx, qpp, 'r', z, pp(ii,:)'*1e-6, 'b');
% xlabel('Position (nm)'); ylabel('Hole concentration (cm^-3)'); 
% title('Hole concentration profile'); grid on;
% % legend set! 
% legend_h = legend('SDevice', 'Matlab');
% set(legend_h,'Location','Best','FontSize',11);
% % error plot! 
% figure('Color', [1,1,1]); % z-slice
% plot(qx, abs( (qpp - pp(ii,:)'*1e-6 )./qpp ), '--ro');
% xlabel('Position (nm)'); ylabel('Hole concentration error'); 
% title('Hole concentration error profile'); grid on;
% legend_h = legend('abs((SDevice-Matlab)/SDevice)');
% set(legend_h,'Location','Best','FontSize',11);
% % find maximum error! 
% max(max(abs( (qpp - pp(ii,:)'*1e-6)./qpp )))

%% mob plot 
        %% position grid setting 
        x = csvread('xmesh.csv', 1, 0, [1 0 69 0])*1e+3;   % node position [nm]
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

        z = csvread('zmesh.csv', 1, 0, [1 0 21 0])*1e+3;
        z_nod = size(z,1);
        z_dlt = (z(2:end) - z(1:end-1));
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

% step #0. node space setting (common w/ hole)
delta_x_n = x_dlt(1:end-1);
delta_x_p = x_dlt(2:end);
delta_z_n = z_dlt(1:end-1)';  
delta_z_p = z_dlt(2:end)';

nx = size(phi, 1);
nz = size(phi, 2);
deltaPhi = zeros(nx, nz, 4);
deltaPhi(2:end,:,  1)   = phi(2:end,  :) - phi(1:end-1,:);    % left  : V(i,j) - V(i-1,j)
deltaPhi(1:end-1,:,2)   = phi(1:end-1,:) - phi(2:end,  :);    % right : V(i,j) - V(i+1,j)
deltaPhi(:,2:end,  3)   = phi(:,2:end  ) - phi(:,1:end-1);    % lower : V(i,j) - V(i,j-1)
deltaPhi(:,1:end-1,4)   = phi(:,1:end-1) - phi(:,2:end  );    % upper : V(i,j) - V(i,j+1)

% - interface considering 
deltaPhi(x_int1 + 2,:,1) = phi(x_int1 + 2,:) - phi(x_int1,    :);
deltaPhi(x_int1,    :,2) = phi(x_int1,    :) - phi(x_int1 + 2,:);
deltaPhi(:,z_int1 + 2,3) = phi(:,z_int1 + 2) - phi(:,    z_int1);
deltaPhi(:,z_int1,    4) = phi(:,    z_int1) - phi(:,z_int1 + 2);
deltaPhi(:,z_int2(1),4) = deltaPhi(:,z_int1(1),4);

% $$$ start : E-Field (new approach)
E_x1 = zeros(nx, nz);
E_x2 = zeros(nx, nz);
E_z3 = zeros(nx, nz);
E_z4 = zeros(nx, nz);

E_x1(2:end-1, :) = bsxfun(@rdivide, deltaPhi(2:end-1,:,1), delta_x_n*1e-9);
E_x2(2:end-1, :) = bsxfun(@rdivide, deltaPhi(2:end-1,:,2), delta_x_p*1e-9);
E_z3(:, 2:end-1) = bsxfun(@rdivide, deltaPhi(:,2:end-1,3), delta_z_n*1e-9);
E_z4(:, 2:end-1) = bsxfun(@rdivide, deltaPhi(:,2:end-1,4), delta_z_p*1e-9);

E_x1(:, z_int1   ) = E_x1(:, z_int1   ).*z_dlt(z_int2(1))/( z_dlt(z_int1(1)-1) + z_dlt(z_int2(1)) );
E_x2(:, z_int1   ) = E_x2(:, z_int1   ).*z_dlt(z_int1(2)-1)/( z_dlt(z_int1(2)-1) + z_dlt(z_int2(2)) );
E_z3(:, z_int1(1)) = 0;
E_z4(:, z_int1(2)) = 0;

E_field = sqrt((E_x1 - E_x2).^2 + (E_z3 - E_z4).^2); 

E_field(1,  :) = deltaPhi(1,  :,2)/(delta_x_n(1  )*1e-9);
E_field(end,:) = deltaPhi(end,:,1)/(delta_x_p(end)*1e-9);

E_field(:, z_int2(1)) = E_field(:, z_int1(1));  % for interface 

mob = canali_n(abs(E_field))*sqrt(2);

mob(:, z_idx{1}) = 0;  
mob(:, z_idx{3}) = 0;  
% $$$ end : E-Field (new approach)

qx = csv_mob(:,1)*1e+3;
qmob = csv_mob(:,2)*1e-4;
% notice: '??' should be changed! 
figure('Color', [1,1,1]); % z-slice
plot(qx, qmob, 'r', z, mob(ii,:)', 'b');
xlabel('Position (nm)'); ylabel('emobility (m^2/V-s)'); 
title('emobility profile'); grid on;
% legend set! 
legend_h = legend('SDevice', 'Matlab');
set(legend_h,'Location','Best','FontSize',11);
% error plot! 
figure('Color', [1,1,1]); % z-slice
plot(qx, abs( (qmob - mob(ii,:)' )./qmob ), '--ro');
xlabel('Position (nm)'); ylabel('emobility error (m^2/V-s)'); 
title('emobility error profile'); grid on;
legend_h = legend('abs((SDevice-Matlab)/SDevice)');
set(legend_h,'Location','Best','FontSize',11);
% find maximum error! 
max(max(abs( (qmob - mob(ii,:)')./qmob )))

%% error check (using node interpolation) 
% % csv value rearrange! 
% ii = 10;
% csv = Vgs0;
% csv(:,1) = csv(:,1)*1e+3;
% qx = csv(:,1);
% qphi = csv(:,2);
% % notice: '??' should be changed! 
% newphi = interp1(x,phi(:,ii),qx);
% newphi(1) = phi(1,ii);
% newphi(end)= phi(end,ii);
% % notice: '??' should be changed! 
% figure('Color', [1,1,1]); % z-slice
% plot(qx, qphi, 'r', x, phi(:,ii), 'b');
% xlabel('Distance (nm)'); ylabel('Potential (V)'); 
% title('Electric potential slice(z-direction) profile'); grid on;
% % legend set! 
% legend_h = legend('SDevice', 'Matlab');
% set(legend_h,'Location','Best','FontSize',11);
% % error plot! 
% figure('Color', [1,1,1]); % z-slice
% plot(qx, abs(newphi - qphi ), '--ro');
% xlabel('Position (nm)'); ylabel('Potential difference (V)'); 
% title('Electric potential difference(z-direction)'); grid on;
% legend_h = legend('abs(SDevice-Matlab)');
% set(legend_h,'Location','Best','FontSize',11);
% % find maximum error! 
% max(abs(newphi - qphi ))