function [ jbase ] = configueJbase()
%% 함수설명 : Jacobian 행렬을 구성합니다. 
% Jacobian (charge + potential) 중 potential 부분을 구성 

%% 글로벌 변수로부터 상수를 불러옵니다. 
global xmesh;   % x방향 mesh 
global zmesh;   % z방향 mesh 
global const_p; % 위치에 따라 달라지는 상수

% x방향 인덱스, 간격 가져오기 
nx = xmesh.nx;          % 총 노드 
x_int1 = xmesh.int1;    % 인덱스 (계면 왼쪽)
x_int2 = xmesh.int2;    % 인덱스 (계면 오른쪽)
x_dlt  = xmesh.dlt;     % 노드 간격

% z방향 인덱스, 간격 가져오기 
nz = zmesh.nz;          % 총 노드 
z_int1 = zmesh.int1;    % 인덱스 (계면 왼쪽)
z_int2 = zmesh.int2;    % 인덱스 (계면 오른쪽)
z_dlt  = zmesh.dlt;     % 노드 간격

% 상수 불러오기 
eps_r    = const_p.eps_r;
boundary = const_p.boundary;
 
%% 계면 위치에서 상대 유전률 수정 
% z방향 계면 위치에서 왼쪽 방향 유전률 수정   
eps_x = eps_r;
eps_x(:,z_int1) = bsxfun(@times, eps_r(:,z_int1), (z_dlt(z_int1-1)./(z_dlt(z_int1-1)+z_dlt(z_int2)))' ) ...
                + bsxfun(@times, eps_r(:,z_int2), (z_dlt(z_int2)  ./(z_dlt(z_int1-1)+z_dlt(z_int2)))' );
            
eps_left        = eps_x(2:end,:);
eps_left(end,:) = eps_x(end,:)*2;   % Neumann 경계조건 

% z방향 계면 위치에서 오른쪽 방향 유전률 수정    
eps_right           = eps_x(1:end-1,:);
eps_right(x_int1,:) = eps_x(x_int2,:);  % x_int1에서 오른쪽 field는 x_int2의 오른쪽 field 참조 
eps_right(1,:)      = eps_x(1,:)*2; % Neumann 경계조건 

% x방향 계면 위치에서 아래쪽 방향 유전률 수정   
eps_z = eps_r;
eps_z(x_int1,:) = bsxfun(@times, eps_r(x_int1,:), x_dlt(x_int1-1)./(x_dlt(x_int1-1)+x_dlt(x_int2)) ) ...
                + bsxfun(@times, eps_r(x_int2,:), x_dlt(x_int2)  ./(x_dlt(x_int1-1)+x_dlt(x_int2)) );
            
eps_lower        = eps_z(:,2:end);
eps_lower(:,end) = eps_z(:,end)*2;  % Neumann 경계조건 

% x방향 계면 위치에서 위쪽 방향 유전률 수정   
eps_upper           = eps_z(:,1:end-1);
eps_upper(:,z_int1) = eps_z(:,z_int2);  % z_int1에서 위쪽 field는 z_int2의 위쪽 field 참조 
eps_upper(:,1)      = eps_z(:,1)*2; % Neumann 경계조건 

%% 각 위치에서 왼쪽, 오른쪽, 아래쪽, 위쪽 field 계산 
% 설명: Poisson 방정식을 적분식으로 변경, 3D로부터 homogeneous방향을 고려한 2D 변환. 
%      즉, 6방향 면적분이 4방향 면적분으로 변환됨 
% 왼쪽 방향 면적분 구성 
delta_n = x_dlt; 
delta_p = [x_dlt(2:end) ; x_dlt(end)];
j_left  = zeros(nx, nz);
j_left(2:end,:)  = bsxfun(@rdivide, eps_left, delta_n.*(delta_p+delta_n)/2 );
% j_left(2:end,:)  = bsxfun(@rdivide, eps_left, 1);
j_cent1 = -j_left;

% 오른쪽 방향 면적분 구성 
delta_n = [x_dlt(1) ; x_dlt(1:end-1)]; 
delta_p = x_dlt;
j_right  = zeros(nx, nz);
j_right(1:end-1,:)  = bsxfun(@rdivide, eps_right, delta_p.*(delta_p+delta_n)/2 );
% j_right(1:end-1,:)  = bsxfun(@rdivide, eps_right, 1);
j_cent2 = -j_right;

% 아래쪽 방향 면적분 구성 
delta_n = z_dlt;
delta_p = [z_dlt(2:end) ; z_dlt(end)];
j_lower  = zeros(nx, nz);
j_lower(:,2:end) = bsxfun(@rdivide, eps_lower, (delta_n.*(delta_p+delta_n)/2)' );
% j_lower(:,2:end) = bsxfun(@rdivide, eps_lower, 1);
j_cent3 = -j_lower;

% 위쪽 방향 면적분 구성 
delta_n = [z_dlt(1) ; z_dlt(1:end-1)];
delta_p = z_dlt;
j_upper  = zeros(nx, nz);
j_upper(:,1:end-1) = bsxfun(@rdivide, eps_upper, (delta_p.*(delta_p+delta_n)/2)' );
% j_upper(:,1:end-1) = bsxfun(@rdivide, eps_upper, 1);
j_cent4 = -j_upper;

%% 계산된 field를 토대로 Jacobian 행렬 구성 
% Jacobian 행렬 저장공간 생성 
jbase = zeros(nx*nz, nx*nz, 2); % 3번째 index는 1: x방향 field, 2: z방향 field를 의미

% x방향으로 행렬 구성 
j_cent1 = matrixToVector(j_cent1,nx,nz);  % 중앙1 : 2차원에서 1차원으로 
j_left  = matrixToVector(j_left, nx,nz);  % 중앙1 왼쪽 : 2차원에서 1차원으로 
j_cent2 = matrixToVector(j_cent2,nx,nz);  % 중앙2 : 2차원에서 1차원으로 
j_right = matrixToVector(j_right,nx,nz);  % 중앙2 오른쪽 : 2차원에서 1차원으로 

% 해당 위치에 field 성분 더함 
jbase(2:end, 2:end,     1) = jbase(2:end, 2:end,     1) + diag(j_cent1(2:end));     % 중앙1(digaonal) 위치
jbase(2:end, 1:end-1,   1) = jbase(2:end, 1:end-1,   1) + diag(j_left (2:end));     % 중앙1 왼쪽 위치
jbase(1:end-1, 1:end-1, 1) = jbase(1:end-1, 1:end-1, 1) + diag(j_cent2(1:end-1));   % 중앙2(digaonal) 위치
jbase(1:end-1, 2:end,   1) = jbase(1:end-1, 2:end,   1) + diag(j_right(1:end-1));   % 중앙2 오른쪽 위치

% z방향으로 행렬 구성 
j_cent3 = matrixToVector(j_cent3,nx,nz);  % 중앙3 : 2차원에서 1차원으로 
j_lower = matrixToVector(j_lower,nx,nz);  % 중앙3 아래쪽 : 2차원에서 1차원으로 
j_cent4 = matrixToVector(j_cent4,nx,nz);  % 중앙4: 2차원에서 1차원으로 
j_upper = matrixToVector(j_upper,nx,nz);  % 중앙4 윈족 : 2차원에서 1차원으로 

% 해당 위치에 field 성분 더함 
jbase(1+nx:end, 1+nx:end, 2) = jbase(1+nx:end, 1+nx:end, 2) + diag(j_cent3(1+nx:end));  % 중앙3(digaonal) 위치
jbase(1+nx:end, 1:end-nx, 2) = jbase(1+nx:end, 1:end-nx, 2) + diag(j_lower(1+nx:end));  % 중앙3 아래쪽 위치
jbase(1:end-nx, 1:end-nx, 2) = jbase(1:end-nx, 1:end-nx, 2) + diag(j_cent4(1:end-nx));  % 중앙4(digaonal) 위치
jbase(1:end-nx, 1+nx:end, 2) = jbase(1:end-nx, 1+nx:end, 2) + diag(j_upper(1:end-nx));  % 중앙4 아래쪽 위치 

% Jacobian 행렬 구성 
jbase = jbase(:,:,1) + jbase(:,:,2);

%% 예외 처리 
% #1. 계면 처리
% 설명 : 두 물질이 만나는 계면이 있다고 하면 계면에 아주 가까운 좌우(또는 상하) 두 개의 노드를 할당하여 해석. 
%        그런데 Poisson 방정식 해석 시에는 두 개가 아닌 하나의 노드에서만 해석 
%        그래서 x방향(x에 수직하는) 계면에는 왼쪽 노드를, z방향(z에 수직하는) 계면에서는 아래쪽 노드를 선택하여 해석. 
%        이에 따라 발생하는 예외 상황들을 직접 수정하였음. 
% 계면 처리 (x방향) : x_int1에서 field를 계산할 때 계면 노드가 아닌 그 오른쪽 노드를 참조하도록 수정. x_int2가 x_int1 복사.   
for i = 1:size(x_int1,1)
    for j = 0:nz-1
        x1 = x_int1(i) + j*nx;
        x2 = x_int2(i) + j*nx;
        jbase(x1,x2+1) = jbase(x1,x2);  % x_int1에서의 방정식에 x_int2의 오른쪽 노드를 참조
        jbase(x1,x2  ) = 0;             % x_int1에서의 방정식에 x_int2를 참조하지 않음 

        jbase(x2,:)  = 0;   % x_int2에서 기존의 식들을 모두 지워버림 
        jbase(x2,x2) = 1;   % x_int2에서의 방정식을 x_int1을 가져오는 식으로 대체 
        jbase(x2,x1) = -1;
    end 
end

% 계면 처리 (z방향) : z_int1에서 field를 계산할 때 계면 노드가 아닌 그 위쪽 노드를 참조하도록 수정. z_int2가 z_int1 복사.
for i = 1:nx
    for j = 1:size(z_int1,1)
        z1 = i+(z_int1(j)-1)*nx;
        z2 = i+(z_int2(j)-1)*nx;
        jbase(z1,z2+nx) = jbase(z1,z2);     % z_int1에서의 방정식에 z_int2의 위쪽 노드를 참조
        jbase(z1,z2   ) = 0;                % z_int1에서의 방정식에 z_int2을 참조하지 않음 

        jbase(z2,:)  = 0;   % z_int2에서 기존의 식들을 모두 지워버림 
        jbase(z2,z2) = 1;   % z_int2에서의 방정식을 z_int1을 가져오는 식으로 대체 
        jbase(z2,z1) = -1;
    end 
end

% 계면 처리 (x방향, z방향 교차 지점) : 대각선에 위치한 xz_int2가 xz_int1 복사
for i = 1:size(x_int1,1)
   for j  = 1:size(z_int1,1)
       xz1 = x_int1(i) + (z_int1(j)-1)*nx;
       xz2 = x_int2(i) + (z_int2(j)-1)*nx;
       
       jbase(xz2, :) = 0;       % xz_int2에서 기존의 식들을 모두 지워버림 
       jbase(xz2, xz2) = 1;     % xz_int2에서 방정식을 xz_int1을 가져오는 식으로 대체 
       jbase(xz2, xz1) = -1;
   end
end

% #2. 경게 조건 처리 
% 경계 조건 설정 (Dirichlet boundary) 
jbase_BC = eye(nx*nz);          
BC_index = find(boundary ~= 0); % 경계 조건에 해당하는 위치를 계산
jbase(BC_index,:) = jbase_BC(BC_index,:);    % identity 행렬에서 해당하는 행을 복사하여 Jacobian 행에 붙여넣기

end






















