function [phi, nn] = initPoisson2D(iterNum, jbase, phi)
%% 함수설명 : 초기 포텐셜 guess를 위한 Poisson 방정식을 해석합니다.
% 각 파라미터는 다음과 같습니다. 
% iterNum : 최대 iteration 횟수
% jbase : Jacobian 행렬 (potential에 대한 정보만 포함)
% phi : 포텐셜

%% 글로벌 변수로부터 상수를 불러옴
global xmesh; % x방향 mesh 
global zmesh; % z방향 mesh 
global const_i; % 위치와 무관한 상수
global const_p; % 위치에 따라 달라지는 상수 

% x방향 인덱스, 간격 가져오기 
x_int1 = xmesh.int1;    % 인덱스 (경계 왼쪽)
x_int2 = xmesh.int2;    % 인덱스 (경계 오른쪽)
x_dlt = xmesh.dlt;      % 노드 간격

% z방향 인덱스, 간격 가져오기 
z_int1 = zmesh.int1;    % 인덱스 (경계 왼쪽)
z_int2 = zmesh.int2;    % 인덱스 (경계 왼쪽)
z_dlt = zmesh.dlt;      % 노드 간격

% 상수 불러오기 
q = const_i.q;
Vt = const_i.Vt;
eps0 = const_i.eps0;
ni = const_p.ni;
doping = const_p.doping;        % 도핑정보 
boundary = const_p.boundary;    % 경계조건 정보
BC_index = find(boundary ~= 0); % 경계조건에 해당하는 위치 index

%% Newton-Raphson iteration을 시작
for i = 1:iterNum
    %% #1. residual vector 계산 
    % charge term을 계산 ρ = q(p - n + Nd+)
    nn = ni.*exp(+phi/Vt);      % 전자농도
    pp = ni.*exp(-phi/Vt);      % 홀농도 
    g  = q*(pp-nn+doping);      % charge 계산 
    
    % 서로 다른 두 물질 계면에서의 예외처리 (계면에서는 두 노드가 할당되어 있음) :
    % <x축과 수직한 계면> (계면에 아주 가까운 왼쪽/오른쪽 두 노드가 할당되어있음)
    % 왼쪽 노드: 다른 물질 계면에서 charge를 계산할 때 계면 왼쪽/오른쪽 노드의 charge를 노드 간격에 따라 가중합
    % 오른쪽 노드: 0 (오른쪽 노드에서는 Poisson 방정식을 해석하지 않음. charge는 반드시 0)
    % <z축과 수직한 계면> (계면에 아주 가까운 아래쪽/위쪽 두 노드가 할당되어있음)
    % 아래쪽 노드: 다른 물질 계면에서 charge를 계산할 때 계면 아래쪽/위쪽 노드의 charge를 노드 간격에 따라 가중합
    % 위쪽 노드: 0 (위쪽 노드에서는 Poisson 방정식을 해석하지 않음. charge는 반드시 0)
    g(x_int1,:) = bsxfun( @times, g(x_int1,:), x_dlt(x_int1-1)./(x_dlt(x_int1-1)+x_dlt(x_int1)) ) ...
                + bsxfun( @times, g(x_int2,:), x_dlt(x_int1)  ./(x_dlt(x_int1-1)+x_dlt(x_int1)) );
    g(:,z_int1) = bsxfun( @times, g(:,z_int1), (z_dlt(z_int1-1)./(z_dlt(z_int1-1)+z_dlt(z_int1)))' ) ...
                + bsxfun( @times, g(:,z_int2), (z_dlt(z_int1)  ./(z_dlt(z_int1-1)+z_dlt(z_int1)))' );
    g(x_int2,:) = 0;    
    g(:,z_int2) = 0;
    
    % Dirichlet 경계 예외처리 (Poisson 방정식을 해석하지 않음. charge는 반드시 0)
    g(BC_index) = 0;

    % 2차원으로 존재하는 행렬을 1차원 Residual vector로 변환 
    % matrixToVector(행렬, x노드갯수, z노드갯수) : 2차원 행렬을 1차원 어레이로 변환하는 함수 
    r = ( jbase*matrixToVector(phi, xmesh.nx, zmesh.nz) - matrixToVector(boundary, xmesh.nx, zmesh.nz) )...
        *eps0*1e+9^2;              % 단위를 고려. (유전률, nm^-2)
    R = r + matrixToVector(g, xmesh.nx, zmesh.nz);

    %% #2. Jacobian matrix 구성 
    % potential에 대해 미분된 charge term을 계산
    h = -(q/Vt)*(pp+nn);
    
    % 다른 물질 계면에서의 예외처리를 #1에서와 동일한 방법으로 처리 
    h(x_int1,:) = bsxfun( @times, h(x_int1,:), x_dlt(x_int1-1)./(x_dlt(x_int1-1)+x_dlt(x_int1)) ) ...
                + bsxfun( @times, h(x_int2,:), x_dlt(x_int1)  ./(x_dlt(x_int1-1)+x_dlt(x_int1)) );
    h(:,z_int1) = bsxfun( @times, h(:,z_int1), (z_dlt(z_int1-1)./(z_dlt(z_int1-1)+z_dlt(z_int1)))' ) ...
                + bsxfun( @times, h(:,z_int2), (z_dlt(z_int1)  ./(z_dlt(z_int1-1)+z_dlt(z_int1)))' );
    h(x_int2,:) = 0;
    h(:,z_int2) = 0;
    
    % Dirichlet 경계 예외처리
    h(BC_index) = 0;
    
    j = jbase*(eps0)*1e+9^2;        % 단위를 고려. (유전률, nm^-2)
    J = j + diag( matrixToVector(h, xmesh.nx, zmesh.nz) );
    
    %% #3. 포텐셜 update 값을 계산 및 수렴성 체크 
    % Jacobian*delta_phi = -Residual -> delta_phi =  -inverse(Jacobian)*Residual을 계산
    invJ = inv(J);
    dphi = -invJ*R;
    dphiMat = vectorToMatrix(dphi, xmesh.nx, zmesh.nz); % 계산된 delta_phi를 1차원 어레이에서 2차원 행렬로 변환
    
    phi = phi + dphiMat;    % 업데이트를 반영 
    dphiVec = matrixToVector(dphiMat, xmesh.nx, zmesh.nz);
    stop(i) = full(max(abs(dphiVec)));  % update 값 중 최대값을 error로 간주하여 확인 
    
%     disp(sprintf('initial Poisson trial[%d]-error: %d \n', i, stop(i)))

    if stop(i) < 1e-12 % error가 기준값 1e-12보다 작으면 수렴판정하여 계산을 종료 
        break;
    end
end

end 
