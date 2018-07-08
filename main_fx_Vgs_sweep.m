% main 스크립트와 구성이 동일합니다. 설명을 위해서는 main.m 을 참조해주세요. 
function [  ] = main_fx_Vgs_sweep( Vds_bias )

% clear;
% set(0,'defaultfigurecolor',[1 1 1]);

%% sequence 
% #0. Define basic information (mesh size / position / constant) 
% #1. Solve 2D init Poisson eq.
% #2. Repeat following loop
%   #2-1. Solve Poisson eq, Current eq (electron, hole)
%   #2-2. If error is larger than 1e-5, repeat #2. 
% #3. Print result such as concentration, potential

% do_init = 0;

% if do_init == 1
    clear -global;
    %% model parameters 
    Vg_bias = 0.0; 
%     Vgate  = 0.533744;###
    Vgate  = 0.5050;
    Vgs = Vgate + Vg_bias; 
%     Vsource = 0.6113137; ###
%     Vss = Vsource;###
%     Vds = Vsource;###
    
    Nd = [2e+20 ; 0 ; 2e+20]; 
    % load mesh information from (.csv) files
    mesh_Generation();

    %% constant setting 
    [ phi ] = initConstSet( Nd, Vgs );
                    
    %% initial Poisson eq. 
    jbase = configueJbase();
%     return;

    [ phi, nn  ] = initPoisson2D( 100, jbase, phi);
    
    load(sprintf('Gate\\result_%03d.mat',Vds_bias*100+1));
%     load(sprintf('Gate\\Vds_%.1f\\result_%03d.mat',Vds_bias, 152));

% ### variable explorer code ###
%     phi = originVariable(1,1,phi);
%     nn = originVariable(1,1,nn);
%     global const_p
%     doping = const_p.doping;
%     doping = originVariable(1,1,doping);
    
%     global xmesh
%     x = xmesh.node;
%     global zmesh
%     z = zmesh.node';
    
    
%     return;
% return;
%     save('init.mat');
%     return;
% else 
%     load('init.mat');
% end

% bias = 0.5;
% nVds = 51;   
% Vds = linspace(0, bias, nVds);     % simulation target 

nVgs = 51;   
Vgs_delta = linspace(0, -1, nVgs);     % simulation target 
Ids = Vgs_delta*0;
Ids_tn = Ids;

Vds = Vds_bias;

% bias = 0.6;
% nVgs = 4;   
% Vgs = linspace(0, bias, nVgs);     % simulation target 

for i = 1:nVgs
    
    deltaE = 0.5e-3;    % 1 meV
    % delta1, nodenum2, valley

    % valley #1(l,t,t) #2(t,l,t) #3(t,t,l)
    [Em_t, Vm_t] = mode_Confinement( 1, phi);   % for #1, #2 
    [Em_l, Vm_l] = mode_Confinement( 3, phi);   % for #3 

    E1 = min(min(Em_t(:,1)),min(Em_l(:,1)))-0.3;    
    E2 = max(max(Em_t(:,2)),max(Em_l(:,5)))+0.3;
    E1 = (round(E1/deltaE))*deltaE;
    E2 = (round(E2/deltaE))*deltaE;
    totalE = E1:deltaE:E2;

    nodeNum = 1000;
    FF1_t = configue_Ffunction( 1, nodeNum, 0, totalE);  % for #1, #3
    FF2_t = configue_Ffunction( 1, nodeNum, Vds, totalE);    
    FF1_l = configue_Ffunction( 2, nodeNum, 0, totalE);     % for #2 
    FF2_l = configue_Ffunction( 2, nodeNum, Vds, totalE);     % for #2 
    
    for j = 1:100
        phi_old = phi;
        tic;
       %% calculation of electron concentration 
        % valley #1(l,t,t) #2(t,l,t) #3(t,t,l)
        [Em_t, Vm_t] = mode_Confinement( 1, phi);   % for #1, #2 
        [Em_l, Vm_l] = mode_Confinement( 3, phi);   % for #3 
        
        nn1 = phi*0; nn2 = nn1; nn3 = nn1; 
        for k = 1:2
            % valley #1(l,t,t)
            [E_v1, FF1_v1, FF2_v1] = find_Ffunction(Em_t, k, totalE, FF1_t, FF2_t);
            nn1 = nn1 + 2*negf_Transport(1, Em_t, Vm_t, k, E_v1, FF1_v1, FF2_v1);
            % valley #2(t,l,t)
            [E_v2, FF1_v2, FF2_v2] = find_Ffunction(Em_t, k, totalE, FF1_l, FF2_l);
            nn2 = nn2 + 2*negf_Transport(2, Em_t, Vm_t, k, E_v2, FF1_v2, FF2_v2);            
        end
        for k = 1:5
            % valley #3(t,t,l)
            [E_v3, FF1_v3, FF2_v3] = find_Ffunction(Em_l, k, totalE, FF1_t, FF2_t);
            nn3 = nn3 + 2*negf_Transport(3, Em_l, Vm_l, k, E_v3, FF1_v3, FF2_v3);
        end
        
        nn = nn1 + nn2 + nn3;
        
%         save('result.mat');
%         elapsedTime = toc
%         return; 

        [ phi, nn ] = nLinPoisson2D( 100, jbase, phi, nn, Vgs_delta(i));

        elapsedTime(i,j) = toc;

%         disp(sprintf('Elapsed Time(sec): %03d \n', elapsedTime(i,j)));
        
        stop(i,j) = max(max(abs(phi - phi_old)));
        disp(sprintf('[%d]self-consist loop[%d]-error: %d \n', i, j, stop(i,j)));
        if (stop(i,j) < 5e-4) || (j == 100)
            break;
        end
    end 
    %% calculation of current 
    for k = 1:2
        % valley #1(l,t,t)
        [E_v1, FF1_v1, FF2_v1] = find_Ffunction(Em_t, k, totalE, FF1_t, FF2_t);
        [Ids_temp, Ids_tn_temp] = negf_Current(1, Em_t, k, E_v1, FF1_v1, FF2_v1);
        Ids(i) = Ids(i) + 2*Ids_temp;
        Ids_tn(i) = Ids_tn(i) + 2*Ids_tn_temp;
        % valley #2(t,l,t)        
        [E_v2, FF1_v2, FF2_v2] = find_Ffunction(Em_t, k, totalE, FF1_l, FF2_l);
        [Ids_temp, Ids_tn_temp] = negf_Current(2, Em_t, k, E_v2, FF1_v2, FF2_v2);
        Ids(i) = Ids(i) + 2*Ids_temp;
        Ids_tn(i) = Ids_tn(i) + 2*Ids_tn_temp;
    end
    for k = 1:5
        % valley #3(t,t,l)
        [E_v3, FF1_v3, FF2_v3] = find_Ffunction(Em_l, k, totalE, FF1_t, FF2_t);
        [Ids_temp, Ids_tn_temp] = negf_Current(3, Em_l, k, E_v3, FF1_v3, FF2_v3);
        Ids(i) = Ids(i) + 2*Ids_temp;
        Ids_tn(i) = Ids_tn(i) + 2*Ids_tn_temp;
    end
    
    %% save and extract detailed variable  
    if (Vgs_delta(i) - floor(Vgs_delta(i)*10)/10 < 1e-10) %% 0 V, 0.1 V, 0.2 V, 0.3 V, 0.4 V, 0.5 V
        % [ A1, A2, T ] = negf_ShowVariables( valleyNum, Em, k_count, E )
%         % valley #1(l,t,t)
%         [ v1_1_A1, v1_1_A2, v1_1_T, v1_1_Ids ] = negf_ShowVariables( 1, Em_t, 1, totalE, FF1_t, FF2_t );
%         [ v1_2_A1, v1_2_A2, v1_2_T, v1_2_Ids ] = negf_ShowVariables( 1, Em_t, 2, totalE, FF1_t, FF2_t );
%         % valley #2(t,l,t)
%         [ v2_1_A1, v2_1_A2, v2_1_T, v2_1_Ids ] = negf_ShowVariables( 2, Em_t, 1, totalE, FF1_l, FF2_l );
%         [ v2_2_A1, v2_2_A2, v2_2_T, v2_2_Ids ] = negf_ShowVariables( 2, Em_t, 2, totalE, FF1_l, FF2_l );
%         % valley #3(t,t,l)
%         [ v3_1_A1, v3_1_A2, v3_1_T, v3_1_Ids ] = negf_ShowVariables( 3, Em_l, 1, totalE, FF1_t, FF2_t );
%         [ v3_2_A1, v3_2_A2, v3_2_T, v3_2_Ids ] = negf_ShowVariables( 3, Em_l, 2, totalE, FF1_t, FF2_t );
%         [ v3_3_A1, v3_3_A2, v3_3_T, v3_3_Ids ] = negf_ShowVariables( 3, Em_l, 3, totalE, FF1_t, FF2_t );
%         [ v3_4_A1, v3_4_A2, v3_4_T, v3_4_Ids ] = negf_ShowVariables( 3, Em_l, 4, totalE, FF1_t, FF2_t );
%         [ v3_5_A1, v3_5_A2, v3_5_T, v3_5_Ids ] = negf_ShowVariables( 3, Em_l, 5, totalE, FF1_t, FF2_t );
%         
        save(sprintf('Gate\\Vds_%.1f\\result_%03d.mat',Vds_bias, 150+2*i));
%         clear v1_1_A1 v1_1_A2 v1_1_T v1_1_Ids
%         clear v1_2_A1 v1_2_A2 v1_2_T v1_2_Ids
%         clear v2_1_A1 v2_1_A2 v2_1_T v2_1_Ids
%         clear v2_2_A1 v2_2_A2 v2_2_T v2_2_Ids
%         clear v3_1_A1 v3_1_A2 v3_1_T v3_1_Ids
%         clear v3_2_A1 v3_2_A2 v3_2_T v3_2_Ids
%         clear v3_3_A1 v3_3_A2 v3_3_T v3_3_Ids
%         clear v3_4_A1 v3_4_A2 v3_4_T v3_4_Ids
%         clear v3_5_A1 v3_5_A2 v3_5_T v3_5_Ids
    else
        save(sprintf('Gate\\Vds_%.1f\\result_%03d.mat',Vds_bias, 150+2*i));
    end
end 
    
% save('result.mat');



end

