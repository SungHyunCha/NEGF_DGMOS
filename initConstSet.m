function [ phi ] = initConstSet( Nd, Vgs )
%% 함수설명 : 시뮬레이션에 필요한 상수들을 글로벌 변수 형태로 저장합니다. 
% 각 파라미터는 다음과 같습니다. 
% Nd : 도핑 농도 (1x3 어레이에 각 지역의 도핑농도를 저장)
% Vgs : Gate단의 Dirichlet 경계값 

%% 글로벌 변수로부터 상수를 불러옵니다. 
global xmesh;   % x방향 mesh 
global zmesh;   % z방향 mesh 

%% 글로벌 변수에 상수를 저장합니다. 
% ni0 백업
% ni0 = 1.0750038e+10;    % [cm^-3]
% ni0 = 1.0758720e+10;    % [cm^-3]

% 위치와 무관한 상수
global const_i;       % 글로벌 변수 생성           
const_i.q = 1.602192e-19;         % [J/eV] or [C]
const_i.Vt = 1.380662e-23*300/const_i.q;    % [V]
const_i.hBar = 1.054571726e-34;   % [J-s] 
const_i.m0 = 9.10938356e-31;      % [kg] 전자 질량
const_i.eps0 = 8.8542e-12;      % [F/m] 
const_i.ni0 = 1.0758720e+10;    % [cm^-3]
const_i.eps_r_si = 11.7;        % [1]
const_i.eps_r_ox = 3.9;         % [1]
const_i.Egap = 1.11;            % [eV] 실리콘 bandgap

% 유효 전자 질량
global mass;        % 글로벌 변수 생성
ml = 0.98*const_i.m0;   % ml = 0.98m0
mt = 0.19*const_i.m0;   % mt = 0.19m0  
mass.m_x = [ml mt mt];  % x방향 유효 질량 [valley #1, #2, #3]
mass.m_y = [mt ml mt];  % y방향 유효 질량 [valley #1, #2, #3]
mass.m_z = [mt mt ml];  % z방향 유효 질량 [valley #1, #2, #3]

% 위치에 따라 달라지는 상수 
global const_p;  
const_p.eps_r      = zeros(xmesh.nx, zmesh.nz);   % 상대유전률 (left, right, bottom, top)
const_p.ni         = zeros(xmesh.nx, zmesh.nz);     % intrinsic 전자 농도 [m^-3]
const_p.doping     = zeros(xmesh.nx, zmesh.nz);     % 도핑 농도 [m^-3]
const_p.boundary   = zeros(xmesh.nx, zmesh.nz);   % 경계 조건 (phi, n, p)

% 포텐셜 변수 생성 
phi = zeros(xmesh.nx, zmesh.nz); 

% 금속 라인 인덱스 생성 
% source 금속 인덱스
x_source = 1;   
z_source = [zmesh.int1(1) zmesh.idx{2} zmesh.int2(2)];  

% drain 금속 인덱스
x_drain  = xmesh.nx;  
z_drain  = z_source;              

% gate1 금속 인덱스
x_gate1  = [xmesh.int1(1) xmesh.idx{2} xmesh.int2(2)];  
z_gate1  = 1;   

% gate2 금속 인덱스
x_gate2  = x_gate1;                         
z_gate2  = zmesh.nz;  

%% 위치에 따라 달라지는 상수 정의
% 상수 setting : 상대유전률
const_p.eps_r(:,zmesh.idx{1}) = const_i.eps_r_ox;    % [1] in oxide.
const_p.eps_r(:,zmesh.idx{2}) = const_i.eps_r_si;    % [1] in semi.
const_p.eps_r(:,zmesh.idx{3}) = const_i.eps_r_ox;    % [1] in oxide.

% 상수 setting : intrinsic 전자 농도 
const_p.ni(:, zmesh.idx{1}) = 0;             % [m^-3] in oxide 
const_p.ni(:, zmesh.idx{2}) = const_i.ni0*1e+6;    % [m^-3] in semi.
const_p.ni(:, zmesh.idx{3}) = 0;             % [m^-3] in oxide

% 상수 setting : 도핑 농도 
const_p.doping(:, zmesh.idx{1}) = 0;                % in oxide 
const_p.doping(xmesh.idx{1}, zmesh.idx{2}) = Nd(1)*1e+6;   % in semi. (n+)
const_p.doping(xmesh.idx{2}, zmesh.idx{2}) = Nd(2)*1e+6;   % in semi. (i)
const_p.doping(xmesh.idx{3}, zmesh.idx{2}) = Nd(3)*1e+6;   % in semi. (n+)
const_p.doping(:, zmesh.idx{3}) = 0;                % in oxide

% Dirichlet 경계값을 입력 (phi)
const_p.boundary(x_gate1, z_gate1) = Vgs;    % gate (lower) 
const_p.boundary(x_gate2, z_gate2) = Vgs;    % gate (upper)
% boundary(x_source,z_source,1) = Vss;    % source (left)
% boundary(x_drain, z_drain, 1) = Vss;    % drain (right)

% 초기 Poisson 방정식에 사용될 포텐셜 guess 
phi = initPotential(phi, const_p.doping, const_p.ni, const_i.Vt);

% Dirichlet 경계조건이 필요한 위치에 포텐셜 덮어씌움
BC_index = find(const_p.boundary ~= 0);     % 경계값이 nonzero인 인덱스를 찾아서
phi(BC_index) = const_p.boundary(BC_index); % 포텐셜에 해당 값을 덮어씌움 
end

