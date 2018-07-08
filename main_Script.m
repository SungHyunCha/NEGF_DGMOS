%% Vds sweep 샘플코드 (채널 크기가 다른 소자 시뮬레이션)
% 디렉토리 확인 
if exist('Drain', 'dir') == 0
    mkdir('Drain');
end
for i = 1:4
    if exist(sprintf('Drain\\%d',i), 'dir') == 0 
        mkdir(sprintf('Drain\\%d',i));
    end
    if exist(sprintf('Drain\\%d\\Vgs_0.0',i), 'dir') == 0 
        mkdir(sprintf('Drain\\%d\\Vgs_0.0',i));
    end
    if exist(sprintf('Drain\\%d\\Vgs_0.5',i), 'dir') == 0 
        mkdir(sprintf('Drain\\%d\\Vgs_0.5',i));
    end
end

movefile('xmesh_0.5.csv','xmesh_backup.csv');   % xmesh_0.5.csv → xmesh_backup.csv로 변경 (기존 파일 백업)
% xmesh_0.5.csv가 x방향 노드 구성파일로 사용 

movefile('xmesh_1.csv','xmesh_0.5.csv');	% xmesh_1.csv → xmesh_0.5.csv로 변경 (로딩)
main_fx_Vds_sweep(0, 1);    % Vgs = 0.0 V, Vds = 0.00, 0.01, ... , 0.50 V 시뮬레이션 수행
main_fx_Vds_sweep(0.5, 1);  % Vgs = 0.5 V, Vds = 0.00, 0.01, ... , 0.50 V 시뮬레이션 수행
movefile('xmesh_0.5.csv','xmesh_1.csv');    % xmesh_0.5.csv → xmesh_1.csv로 변경 (해제)

movefile('xmesh_2.csv','xmesh_0.5.csv');    % xmesh_2.csv → xmesh_0.5.csv로 변경 (로딩)
main_fx_Vds_sweep(0, 2);
main_fx_Vds_sweep(0.5, 2);
movefile('xmesh_0.5.csv','xmesh_2.csv');    % xmesh_0.5.csv → xmesh_2.csv로 변경 (해제)

movefile('xmesh_3.csv','xmesh_0.5.csv');    % xmesh_3.csv → xmesh_0.5.csv로 변경 (로딩)
main_fx_Vds_sweep(0, 3);
main_fx_Vds_sweep(0.5, 3);
movefile('xmesh_0.5.csv','xmesh_3.csv');    % xmesh_0.5.csv → xmesh_3.csv로 변경 (해제)

movefile('xmesh_4.csv','xmesh_0.5.csv');    % xmesh_4.csv → xmesh_0.5.csv로 변경 (로딩)
main_fx_Vds_sweep(0, 4);
main_fx_Vds_sweep(0.5, 4);
movefile('xmesh_0.5.csv','xmesh_4.csv');    % xmesh_0.5.csv → xmesh_4.csv로 변경 (해제)

movefile('xmesh_backup.csv','xmesh_0.5.csv');   % xmesh_0.5.csv → xmesh_backup.csv로 변경 (기존 파일 복구)

%% Vgs sweep 샘플코드
% bias상태에서 미리 계산된 결과가 필요함 
main_fx_Vgs_sweep(0.1);
main_fx_Vgs_sweep(0.2);
main_fx_Vgs_sweep(0.3);
main_fx_Vgs_sweep(0.4);
main_fx_Vgs_sweep(0.5);