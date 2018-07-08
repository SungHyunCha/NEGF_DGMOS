function [  ] = mesh_Generation(  )
%% 함수설명 : 지정된 위치에서 .csv을 읽어 mesh를 글로벌 변수에 저장합니다. 

%% 위치 노드 설정 
% x노드 불러와서 위치, 간격, 경계를 정의 
last = csvread('xmesh_0.5.csv', 0, 0, [0 0 0 0]);   % 노드 총 갯수 
posOfNode = csvread('xmesh_0.5.csv', 1, 0, [1 0 last 0])*1e+3;   % 노드 위치 [nm]
numberOfNode = size(posOfNode,1);                   % 노드 총 갯수 
delta = (posOfNode(2:end) - posOfNode(1:end-1));    % 노드 간격 [nm]
int1 = find(delta == 0);        % 경계 왼쪽 인덱스 
int2 = int1 + 1;                % 경계 오른쪽 인덱스 
delta(int1) = delta(int2);    

% 경계로 구분된 지역 별 인덱스 모음 (셀 형태로 저장)
idx = {1:int1(1)};            % 지역 1
for i = 2:size(int1,1)        % 사이에 위치한 지역 
    idx{end+1} = int2(i-1):int1(i);
end
idx{end+1} = int2(i):numberOfNode;  % 지역 end

% x정보 글로벌 변수에 저장 
global xmesh;
xmesh.pos = posOfNode;      % 위치 (어레이)
xmesh.nx = numberOfNode;    % 총 노드 갯수 
xmesh.dlt = delta;          % 노드 간격 (어레이)
xmesh.int1 = int1;          % 경계 왼쪽 노드 (어레이)
xmesh.int2 = int2;          % 경계 오른쪽 노드 (어레이)
xmesh.idx = idx;            % 지역 별 인덱스 (어레이)

x = xmesh.pos;
x_int = xmesh.int2;
x(x_int) = [];
xmesh.node = x;             % 위치 (경계 중복 노드 제거)

% z노드 불러와서 위치, 간격, 경계를 정의 
posOfNode = csvread('zmesh_0.125.csv', 1, 0, [1 0 59 0])*1e+3;  % 노드 위치 [nm]
numberOfNode = size(posOfNode,1);                               % 노드 총 갯수
delta = (posOfNode(2:end) - posOfNode(1:end-1));  % 노드 간격 [nm]
%     z_dlt = round(z_dlt*1e+3)*1e-3;
int1 = find(delta == 0);        % 경계 왼쪽 인덱스 
int2 = int1 + 1;                % 경계 오른쪽 인덱스 
delta(int1) = delta(int2);

% 경계로 구분된 지역 별 인덱스 모음 (셀 형태로 저장)
idx = {1:int1(1)};          % 지역 1
for i = 2:size(int1,1)      % 사이에 위치한 지역 
    idx{end+1} = int2(i-1):int1(i);
end
idx{end+1} = int2(i):numberOfNode;      % 지역 end

% z정보 글로벌 변수에 저장 
global zmesh;
zmesh.pos = posOfNode;      % 위치 (어레이)
zmesh.nz = numberOfNode;    % 총 노드 갯수 
zmesh.dlt = delta;          % 노드 간격 (어레이)
zmesh.int1 = int1;          % 경계 왼쪽 노드 (어레이)
zmesh.int2 = int2;          % 경계 오른쪽 노드 (어레이)
zmesh.idx = idx;            % 지역 별 인덱스 (어레이)

z_int = [zmesh.int1(1) zmesh.int2(2)];
z = zmesh.pos;
z(z_int) = [];
zmesh.node = z;
zmesh.int3 = z_int;         % 위치 (경계 중복 노드 제거)

end
