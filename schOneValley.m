function [ nn ] = schOneValley( delta1, nodeNum2, valley, phi, Vds, z, x, x_int2, Egap, Vt, q )
% constant setting 
hBar = 1.054571726e-34;        % [J-s]
% phi_F = 0.560983627;


% node setting (z)
% qz = interp1q((1:size(z,1))', z, (1:0.01:size(z,1))')*1e-9;
z = z*1e-9;
z_dlt = z(2:end) - z(1:end-1);
z_dlt = round(z_dlt*1e+12)*1e-12;

nx = size(phi, 1)-2;
nz = size(phi, 2);
Em = zeros(nx, nz-2);
Vm = zeros(nx, nz-2, nz-2);

% node setting (x)
phi(x_int2,:) = []; % delete interface shadow 
x(x_int2) = []; % delete interface shadow 
x = x*1e-9;
x_dlt = x(2:end) - x(1:end-1);
x_dlt = round(x_dlt*1e+12)*1e-12;


% mass setting 
m0 = 9.10938356e-31;       % electron mass [kg] ml = 0.98m0, mt = 0.19m0  sqrt(0.98*0.19)*
ml = 0.98*m0;
mt = 0.19*m0;

if valley == 1 % direction : z 
    m_x = mt;
    m_y = mt;
    m_z = ml;
    k_count = 5;
elseif valley == 2 % direction : y
    m_x = mt;
    m_y = ml;
    m_z = mt;
    k_count = 2;
elseif valley == 3  % direction : x
    m_x = ml;
    m_y = mt;
    m_z = mt;
    k_count = 2;
else 
    disp(sprintf('option : out of range!'));
    return;
end 

for i = 1:nx

    %% input parameters
    Vsch = q*(-phi(i,2:end-1)+Egap/2); % potential energy in J 
%     Vsch = q*(0*phi(i,2:end-1)+Egap/2);
%     Vsch = q*(zeros(nz-2,1));

    %% calculation STEP1: Calculate Hemiltonian
    T = zeros(nz-2, nz-2);
    dlt_n = z_dlt(1:end-1);
    dlt_p = z_dlt(2:end);
    dlt_m = ( dlt_p + dlt_n )/2;
    
    t_left = diag( -(hBar^2/(2*m_z)).*(1./dlt_n).*(1./dlt_m) );
    t_right = diag( -(hBar^2/(2*m_z)).*(1./dlt_p).*(1./dlt_m) );
    
    % left (general)
    T = T - t_left;
    T(2:end, 1:end-1) = T(2:end, 1:end-1) + t_left(2:end, 2:end);
    % right (general)
    T = T - t_right;
    T(1:end-1, 2:end) = T(1:end-1, 2:end) + t_right(1:end-1, 1:end-1);
    % first & last 
    T(1,1) = -2*T(1,2);
    T(end,end) = -2*T(end,end-1);
    
    H = T + diag(Vsch);
    
    % basis(1:9,1:9)
    [vector, value] = eig(H);  %vector: col is eigenvector, value: diagonal element is eigenvalue
    [value, index] = sort(diag(value));
    vector = vector(:, index);

    
    %% calculation STEP2: wavefunction and eigen-energy
    normal = sum( bsxfun(@times, vector.^2, dlt_m),1 );  % integration
    wave_sim = bsxfun(@rdivide, vector, sqrt(normal));   % normalization 
    Ek_sim = (value)';  %[J]
    
    Em(i,:) = Ek_sim;
    Vm(i,:,:) = wave_sim;
end

i = sqrt(-1);
small_i = i*1e-12*q;

sigma1 = zeros(nx,nx);
sigma2 = zeros(nx,nx);

nn_x = zeros(nx,nz-2);

%% calculation STEP1: Calculate Hemiltonian
   
%% first expression
% integ_sim1 = zeros(nx,1);
% for pos = 1:nx 
%     % define E delta 
%     D = (1/(2*pi))*(sqrt(m_x*m_y)/(hBar^2));
%     % define Energy space 
%     E = linspace(Em(pos,1), Em(pos,1) + 0.5*q, 2001);
%     dE =  E(2:end)-E(1:end-1);
%     E = E(1:end-1) + dE/2;
%     %             E = E(1:end-1);
%     F = 1./(1+exp(E/(q*Vt)));
%     %         save schrodinger;
%     integ_sim1(pos) = sum(dE(1).*F*D);
% end

%% NEGF     
t = +(hBar^2/(2*m_x*x_dlt(1)^2));
sbase_first_row = zeros(nx,1); 
sbase_first_row(1) = -2; sbase_first_row(2) = 1;
sbase = sparse(toeplitz(sbase_first_row, sbase_first_row'));
T = -t*sbase;

a = x_dlt(1);     % 노드 배치 간격 
dE = delta1*q; 

for k = 1:k_count 
    E = min(Em(:,k)) -0.25*q:dE:max(Em(:,k)) + 0.5*q;
%     if k == 1
%         if valley == 1 % valley #2 
%         save schrodinger_1;
%         elseif valley == 2 % valley #3 
%         save schrodinger_2;
%         elseif valley == 3  % valley #1 
%         save schrodinger_3;
%         else 
%         disp(sprintf('option : out of range!'));
%         end
%     end
        
%     dE = E(2) - E(1);
%     E = E(1:end-1) + dE/2;
    nE = size(E,2)-1;

    FF1 = Fhalfinv( -(E(1:end-1) + dE/2)/q, nodeNum2, hBar, a, m_y, Vt, q);
    FF2 = Fhalfinv( -Vds -(E(1:end-1) + dE/2)/q, nodeNum2, hBar, a, m_y, Vt, q);
    % NEGF 해석 시작
    A1 = zeros(nx, nE);
    A2 = zeros(nx, nE);
    for j = 1:nE
        E_l = E(j) + dE/2;  % mid energy grid
        ka = acos( 1 - ( E_l - Em(1,k) )/(2*t) );   
        sigma1(1,1) = -t*exp(i*ka);         % self energy calc
        ka = acos( 1 - ( E_l - Em(nx,k) )/(2*t) );
        sigma2(nx, nx) = -t*exp(i*ka);
        G = ( E_l*eye(nx) -T -diag(Em(:,k)) - sigma1 - sigma2 )\eye(nx, nx);

        A1(:,j) = real(diag(i*G*(sigma1 - sigma1')*G'));
        A2(:,j) = real(diag(i*G*(sigma2 - sigma2')*G'));
    end
    % NEGF 해석 완료
    integ_sim2 = ( sum(bsxfun(@times, A1, FF1),2) + sum(bsxfun(@times, A2, FF2),2) )/(2*pi)*dE;
    for pos = 1:nx
%         save tt
        nn_x(pos,:) = nn_x(pos,:) + 2*integ_sim2(pos)*Vm(pos,:,k).^2;
    end
end

if min(min(nn_x)) < 0
    disp('Error: n has negative value');
end

nn_index = [(1:x_int2(1)-1)';...
            (x_int2(1)+1:x_int2(2)-1)';...
            (x_int2(2)+1:nx+2)'];

nn = zeros(nx+2,nz);
nn(nn_index, 2:end-1) = nn_x;
nn(x_int2, :) = nn(x_int2-1, :);
   
% if valley == 1 % valley #2 
%     save schrodinger_1;
% elseif valley == 2 % valley #3 
%     save schrodinger_2;
% elseif valley == 3  % valley #1 
%     save schrodinger_3;
% else 
%     disp(sprintf('option : out of range!'));


return;
   
%% analytic calculation 
% const setting
L = 5e-9;

% Energy : analytic calculation
Ek_ana = (1:Ek_index).^2*(hBar*pi)^2/(2*m_SCH*L^2);

% Wavefunction : analytic calculation
for n = 1:Ek_index
    wave_ana(:,n) = sqrt(2/L)*sin(n*pi*(z(2:end-1)+L/2)/L);
end

% Concentration : analytic calculation
% define E delta 
Ek_ana_start = Ek_ana;
Ek_ana_end = Ek_ana;
Ek_ana_end(:) = Ek_ana(Ek_index);

% integration 
integ_ana = F(Ek_ana_start,Ek_ana_end)*D;
nk_ana = bsxfun(@times, 2*wave_ana.^2, integ_ana);
for i = 2:nx-1
    nn_ana(i-1,2:end-1) = sum(nk_ana(:,1:Ek_index),2)';
end


end
































