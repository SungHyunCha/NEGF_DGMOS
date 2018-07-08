clear;
set(0,'defaultfigurecolor',[1 1 1]);

load(sprintf('result (Vg = 0.0V)\\result_%03d.mat',51));
% return; 

global xmesh;
x = xmesh.pos;
x(xmesh.int2) = [];
global zmesh;
z = zmesh.pos; 

E_l = totalE(1:end-1) + deltaE/2;

%% Em : Subband minimum 
% % valley #1(l,t,t)
% % valley #2(t,l,t)   
% figure; plot(x, Em_t(:,1));
% hold on; plot(x, Em_t(:,2));
% % valley #3(t,t,l)
% figure; plot(x, Em_l(:,1));
% hold on; plot(x, Em_l(:,2));
% hold on; plot(x, Em_l(:,3));
% hold on; plot(x, Em_l(:,4));
% hold on; plot(x, Em_l(:,5));
% % return;

%% FF1(E_l), FF2(E_l) : Fermi function
% % valley #1(l,t,t)
% % valley #3(t,t,l)
% figure; plot(E_l, FF1_t);
% % valley #2(t,l,t)   
% hold on; plot(E_l, FF1_l);
% return;

%% A1(E_l), A2(E_l) : spectral density 
% [X,Y] = meshgrid(E_l,x);
% figure; surface(X,Y,v1_1_A2);
% return;
% hold on; plot('E_l, v2_1_A1');
% hold on; plot('E_l, v3_1_A1');

%% T(E_l) : transmission coefficient 
% figure; plot(E_l, v1_1_T);
% figure; plot(E_l, v2_1_T+v2_2_T);
% figure; plot(E_l, v3_1_T+v3_2_T+v3_3_T+v3_4_T+v3_5_T);

%% Ids(E_l) : current density 
% figure; plot(E_l, v1_1_Ids+v1_2_Ids);
% hold on; plot(E_l, v2_1_Ids+v2_2_Ids);
% hold on; plot(E_l, v3_1_Ids+v3_2_Ids);
% figure; plot(E_l, v1_1_Ids+v2_1_Ids+v3_1_Ids);

%% n(x,z) : electron concentration 
% figure; plot(xmesh.pos, nn(:,30));

%% phi(x,z) : electrostatic potential 
% figure; plot(xmesh.pos, phi(:,30));

%% Ids : current density 
figure; plot(Vds, Ids);


