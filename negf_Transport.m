function [ nn ] = negf_Transport( valleyNum, Em, Vm, k_count, E, FF1, FF2 )
%% 함수설명 : NEGF를 이용하여 전류값을 계산합니다.
% 각 파라미터는 다음과 같습니다. 
% valleyNum : x방향으로의 valley 번호 (#1: m_l, #2 & #3 : m_t)
% Em : mode 에너지(subband minimum)
% Vm : 파동함수 
% k_count : 해석할 mode 번호 
% E : 해석에 사용할 에너지 노드 어레이 
% FF1 : Source단에서의 Fermi 함수값, E index에 대응
% FF2 : Drain단에서의 Fermi 함수값, E index에 대응

%% 글로벌 상수를 불러옵니다. 
if (3 < valleyNum ) || (valleyNum < 1)  % 범위 체크 (1~3)
    disp(sprintf('option : out of range!'));
    return;
end
global mass;    % valley에 따른 전자 유효 질량 
m_x = mass.m_x(valleyNum);

global xmesh;   % x방향 mesh 
global zmesh;   % z방향 mesh
global const_i; % 위치와 무관한 상수

% x방향 간격, 총 노드 갯수 가져오기
x_dlt = xmesh.dlt(1)*1e-9;
nx = xmesh.nx - 2;
% z방향 실리콘 영역 노드 갯수 가져오기 
nz = size(zmesh.idx{2}, 2);

% 상수 불러오기
q = const_i.q;
hBar = const_i.hBar;

% 에너지 단위 변환 ([eV] -> [J])
E = E*q;
Em = Em(:,k_count)*q;

% 해당 mode 번호에 맞는 파동함수 가져오기 
Vm = Vm(:,:,k_count);

%% NEGF 해석을 수행  
% Hamiltonian 행렬 중 T를 구성 : H = T(kinetic energy) + U(potential energy)
t = +(hBar^2/(2*m_x*x_dlt^2));  % 상수 t = hBar^2/(2 m a^2)
sbase_first_row = zeros(nx,1); 
sbase_first_row(1) = -2; sbase_first_row(2) = 1;
sbase = sparse(toeplitz(sbase_first_row, sbase_first_row'));
T = -t*sbase;   % T 완성

% 에너지 노드 재정리 
nE = size(E,2)-1;   % 계산에 사용될 에너지 노드 갯수, 노드 사이의 중점을 구분구적법에 사용
dE = E(2) - E(1);   % 에너지 노드 간격

% 상수 정의
i = sqrt(-1);        % 허수 i 정의 

% 변수 정의 
sigma1 = zeros(nx,nx);  % self 에너지 정의 (source)
sigma2 = zeros(nx,nx);  % self 에너지 정의 (drain)

nn_x = zeros(nx,nz);    % NEGF를 통해 얻는 2차원 전자농도 

A1 = zeros(nx, nE);     % spectral density (Local DOS) - source
A2 = zeros(nx, nE);     % spectral density (Local DOS) - drain 

% NEGF를 통한 전자농도 계산 시작 
for j = 1:nE
    % 계산할 에너지 선택 
    E_l = E(j) + dE/2;  % 구분구적법 계산을 위해 에너지 노드에서 중간값을 선택 
    
    % 셀프 에너지 계산 
    % matlab acos 함수가 [-1,+1]이 [π, 0]에 대응되는 건 맞으나 
    % [-1,+1] 밖의 영역에서는 복소수가 등장하여 원하지 않는 동작을 함 
    % 개선이 필요할 것으로 보임 
    E_ratio1 = 1 - ( E_l - Em(1) )/(2*t);
    ka1 = acos( E_ratio1 +2*(E_ratio1 < -2) +2*(E_ratio1 < -4) );   
    sigma1(1,1) = -t*exp(i*ka1);        
    E_ratio2 = 1 - ( E_l - Em(nx) )/(2*t);
    ka2 = acos( E_ratio2 +2*(E_ratio2 < -2) +2*(E_ratio2 < -4) );
    sigma2(nx, nx) = -t*exp(i*ka2);
    
    % 그린 함수 계산 
    G = ( E_l*eye(nx) -T -diag(Em) - sigma1 - sigma2 )\eye(nx, nx);

    A1(:,j) = real(diag(G*i*(sigma1 - sigma1')*G'));    % source로부터의 Local DOS 
    A2(:,j) = real(diag(G*i*(sigma2 - sigma2')*G'));    % draindm로부터의 Local DOS 
end

% NEGF 해석 후 에너지에 따라 분포하는 전자농도를 더함 (에너지 적분)
integ_sim2 = ( sum(bsxfun(@times, A1, FF1),2) + sum(bsxfun(@times, A2, FF2),2) )/(2*pi)*dE;
for pos = 1:nx % 각 x위치에 대하여 
    % 파동함수의 제곱을 곱하여 z방향으로의 전자농도 확장 
    % factor 2는 전자 스핀
    nn_x(pos,:) = nn_x(pos,:) + 2*integ_sim2(pos)*Vm(pos,:).^2;
end

if min(min(nn_x)) < 0
    disp('Error: n has negative value');
end

% 전자농도에 대하여 x축과 수직 방향의 계면에 노드를 하나씩 추가. (potential과 호환을 위함)
nn_index = [(1:xmesh.int2(1)-1)';...
            (xmesh.int2(1)+1:xmesh.int2(2)-1)';...
            (xmesh.int2(2)+1:nx+2)'];

nn = zeros(nx+2, zmesh.nz);
nn(nn_index, zmesh.idx{2}) = nn_x;
nn(xmesh.int2, :) = nn(xmesh.int2-1, :);

end

