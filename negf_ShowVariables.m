function [ A1, A2, TT, Ids ] = negf_ShowVariables( valleyNum, Em, k_count, E, FF1, FF2 )
%% 함수설명 : NEGF를 이용하여 세부 정보를 계산합니다. 
% 각 파라미터는 다음과 같습니다. 
% valleyNum : x방향으로의 valley 번호 (#3: m_l, #1 & #2 : m_t)
% Em : mode 에너지(subband minimum)
% k_count : 해석할 mode 번호 
% E : 해석에 사용할 에너지 노드 어레이 
% FF1 : Source단에서의 Fermi 함수값, E index에 대응
% FF2 : Drain단에서의 Fermi 함수값, E index에 대응

%% 글로벌 상수를 불러옵니다. 
if (3 < valleyNum ) || (valleyNum < 1) % 범위 체크 (1~3)
    disp(sprintf('option : out of range!'));
    return;
end
global mass;        % valley에 따른 전자 유효 질량 
m_x = mass.m_x(valleyNum);

global xmesh;   % x방향 mesh
global const_i;	% 위치와 무관한 상수

% x방향 노드 간격, 총 노드 갯수 가져오기
x_dlt = xmesh.dlt(1)*1e-9;
nx = xmesh.nx - 2;

% 상수 불러오기
q = const_i.q;
hBar = const_i.hBar;

% 에너지 단위 변환 ([eV] -> [J])
E = E*q;
Em = Em(:,k_count)*q;


%% NEGF 해석을 수행  
% Hamiltonian 행렬 중 T를 구성 : H = T(kinetic energy) + U(potential energy)
t = +(hBar^2/(2*m_x*x_dlt^2));      % 상수 t = hBar^2/(2 m a^2)
sbase_first_row = zeros(nx,1); 
sbase_first_row(1) = -2; sbase_first_row(2) = 1;
sbase = sparse(toeplitz(sbase_first_row, sbase_first_row'));
T = -t*sbase;     % T 완성

% 에너지 노드 재정리 
nE = size(E,2)-1;   % 계산에 사용될 에너지 노드 갯수, 노드 사이의 중점을 구분구적법에 사용
dE = E(2) - E(1);   % 에너지 노드 간격

% 상수 정의
i = sqrt(-1);        % 허수 i 정의 

% 변수 정의 
sigma1 = zeros(nx,nx);  % self 에너지 정의 (source)
sigma2 = zeros(nx,nx);  % self 에너지 정의 (drain)

% 세부 정보를 저장할 변수 
A1 = zeros(nx, nE);     % spectral density (Local DOS) - source
A2 = zeros(nx, nE);     % spectral density (Local DOS) - drain 
TT = zeros(1, nE);      % longitudinal 에너지에 따른 투과계수 T 
Ids = zeros(1, nE);     % longitudinal 에너지에 따른 전류

% NEGF를 통한 전자농도 계산 시작 
for j = 1:nE
    % 계산할 에너지 선택 
    E_l = E(j) + dE/2;      % 구분구적법 계산을 위해 에너지 노드에서 중간값을 선택 
    
    % 셀프 에너지 계산 
    % matlab acos 함수가 [-1,+1]이 [π, 0]에 대응되는 건 맞으나 
    % [-1,+1] 밖의 영역에서는 복소수가 등장하여 원하지 않는 동작을 함 
    % 개선이 필요할 것으로 보임 
    E_ratio1 = 1 - ( E_l - Em(1) )/(2*t);
    ka1 = acos( E_ratio1 +2*(E_ratio1 < -2) +2*(E_ratio1 < -4) );   
    sigma1(1,1) = -t*exp(i*ka1);         % self energy calc
    E_ratio2 = 1 - ( E_l - Em(nx) )/(2*t);
    ka2 = acos( E_ratio2 +2*(E_ratio2 < -2) +2*(E_ratio2 < -4) );
    sigma2(nx, nx) = -t*exp(i*ka2);
    
    % 그린 함수 계산 
    G = ( E_l*eye(nx) -T -diag(Em) - sigma1 - sigma2 )\eye(nx, nx);
    
    A1(:,j) = real(diag(G*i*(sigma1 - sigma1')*G'));    % source로부터의 Local DOS 
    A2(:,j) = real(diag(G*i*(sigma2 - sigma2')*G'));    % draindm로부터의 Local DOS 
    
    gamma1 = i*(sigma1 - sigma1');      % 감마(셀프에너지 허수파트) (source)
    gamma2 = i*(sigma2 - sigma2');      % 감마(셀프에너지 허수파트) (drain)
    TT(j) = sum(real(diag(gamma1*G*gamma2*G')));    % 투과 계수 T
    Ids(j) = 2*q*x_dlt/(2*pi*hBar)*(FF1(j)-FF2(j))*TT(j)*q; % factor 2는 전자 스핀
end

end

