function [ Em, Vm ] = mode_Confinement( valleyNum, phi )
%% 함수설명 : 각각의 x위치에서 z방향으로 Schrodinger 방정식을 해석, mode 에너지와 파동함수를 얻습니다. 
% 각 파라미터는 다음과 같습니다. 
% valleyNum : z방향으로의 valley 번호 (#1 & #2 : m_t, #3: m_l)
% phi : Schrodinger 방정식 해석에 사용할 포텐셜 

%% 글로벌 상수를 불러옵니다. 
if (3 < valleyNum ) || (valleyNum < 1)  % 범위 체크 (1~3)
    disp(sprintf('option : out of range!'));
    return;
end
global mass;    % valley에 따른 전자 유효 질량 
m_z = mass.m_z(valleyNum);

global xmesh;   % x방향 mesh
global zmesh;   % z방향 mesh
global const_i; % 위치와 무관한 상수

% 노드 간격, 총 노드 갯수 가져오기
z_dlt = zmesh.dlt(zmesh.idx{2}(1))*1e-9;
nx = xmesh.nx - 2;
nz = size(zmesh.idx{2}, 2);

% 모드 에너지, 파동함수 저장공간 확보
Em = zeros(nx, nz);
Vm = zeros(nx, nz, nz);

% 노드 변환 : 
% 다음 계산 과정인 NEGF을 고려하여 계면에 2개 할당된 노드 중 한 개를 삭제하고 oxide영역 제거
phi(xmesh.int2,:) = [];     % 계면 오른쪽 노드를 삭제 
phi = phi(:, zmesh.idx{2}); % 실리콘 영역의 포텐셜만 가져옴 

% 상수 불러오기
q = const_i.q;
hBar = const_i.hBar;

%% Schrodinger 방정식 해석을 시작 
for i = 1:nx

    %% 포텐셜을 conduction 밴드 에너지로 변환 
    phi_F = 0.560983627;
    Vsch = q*(-phi(i,:)+phi_F); % potential energy in J 

    % Hamiltonian 행렬 중 T를 구성 : H = T(kinetic energy) + U(potential energy)
    t = +(hBar^2/(2*m_z)).*(1./z_dlt^2);
    sbase_first_row = zeros(nz,1); 
    sbase_first_row(1) = -2; sbase_first_row(2) = 1;
    sbase = sparse(toeplitz(sbase_first_row, sbase_first_row'));
    T = -t*sbase;   % T 완성

    H = T + diag(Vsch); % Hamiltonian 행렬 구성
    
    % Eigen value problem을 해석 및 정렬
    [vector, value] = eig(H);   % vector : Eigen 벡터, value : Eigen 에너지 
    [value, index] = sort(diag(value));     % Eigen 에너지를 오름차순으로 정렬 
    vector = vector(:, index);             

    % mode 에너지와 파동함수를 계산 
    normal = sum( bsxfun(@times, vector.^2, z_dlt),1 );  % normalization factor 계산 
    wave_sim = bsxfun(@rdivide, vector, sqrt(normal));   % normalization을 통해 파동함수 계산
    Ek_sim = (value)';  %[J]
    
    Em(i,:) = Ek_sim/q;     % mode 에너지 [eV]
    Vm(i,:,:) = wave_sim;   % 파동함수 [m^-0.5]
end



end

