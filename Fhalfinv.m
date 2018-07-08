function [ integration ] = Fhalfinv( E, nE, m)
%% 주어진 에너지에 따른 Fermi function을 계산합니다. (y(homogenous) 방향으로 Fermi-Dirac 적분)
% 각 파라미터는 다음과 같습니다. 
% E: 에너지 영역 (1차원 어레이만 처리 가능)
% nE : 에너지 노드 갯수
% m : 유효 전자 질량 (homogeneous 방향)

% 행벡터는 열벡터로 변환하여 처리. 
tOption = 0; % transpose option 
if size(E,2) == 1   % 열벡터 
    tOption = 0;
elseif size(E,1) == 1   % 행벡터
    tOption = 1;
    E = E';
else
    disp('Error: Invalid Input Format (Fermi integral)');
    integration = 0;
    return;
end

% 글로벌 변수로부터 상수 불러오기 
global xmesh;
delta = xmesh.dlt(1)*1e-9;  % nm 환산 
global const_i;
q = const_i.q;
Vt = const_i.Vt;
hBar = const_i.hBar;

% kBT로 나누어진 Ey를 기준으로 적분 
eps = linspace(0, 10, nE);  % y방향 에너지 노드 [1]
d_eps = eps(2) - eps(1);    % 에너지 간격 
eps = eps(1:end-1) + d_eps/2;   % 적분(구분구적법)을 위해 노드 간 중간값을 사용 
[eps_mat, E_mat] = meshgrid(eps, E);    % 2차원으로 전개하여 계산 
% (열이 달라지면 y방향 에너지가, 행이 달라지면 주어진 에너지가 변함)
eps2_mat = eps_mat + d_eps/2;   % 적분 막대기 왼쪽 위치 
eps1_mat = eps_mat - d_eps/2;   % 적분 막대기 오른쪽 위치 

% 에너지 적분을 계산 
integration = (hBar*delta)^-1*sqrt(m*q*Vt/2)/pi*sum( 2*(1./(1+exp(eps_mat - E_mat/Vt))).*(sqrt(eps2_mat) - sqrt(eps1_mat)), 2 );

if tOption == 1 % 입력이 행벡터였으면 transpose
    integration = integration';
end 
end

