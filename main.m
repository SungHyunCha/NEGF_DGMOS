clear;
set(0,'defaultfigurecolor',[1 1 1]);

%% 시뮬레이션 수행 절차 
% #0. 기본 정보를 정의 (메쉬 사이즈 / 위치 / 상수 정의) 
% #1. 2D Poisson 방정식을 해석 (초기 potential 분포를 guess)
% #2. 평형상태에서 다음 self consistent loop를 수행 
%   #2-1. NEGF 해석 (입력: potential, 출력: 전자농도)
%   #2-2. Poisson 방정식 해석 (입력: 전자농도, 출력: potential)
%   #2-3. (#2-1, #2-2 실행 전후 potential 변화 최대값을 error로 정의) 
%         i) error가 0.5e-5 보다 작으면 수렴으로 간주하여 #2-4로 이동
%         ii) error가 0.5e-5 보다 크면 #2-1부터 다시 수행 
%   #2-4. 수렴한 solution(potential)으로부터 전류를 계산, solution을 저장 
% #3. 계산할 다음 bias point를 step만큼 조정하여 #2. 를 수행

%% 초기 변수 설정 
    clear -global;  % 기존 글로벌 변수를 모두 삭제 
    %% bias 정의
    Vg_bias = 0.1;  % Vgs [V]
    % 금속 - workfunction: 4.1 eV 
    % 실리콘 - electron affinity: 4.05 eV, bandgap: 1.11 eV
    % barrier height = 4.1 - (4.05 + 1.11/2) = -0.505 eV -> 0.505 V
    Vg_barrier  = 0.505; % [V]  
    Vgs = Vg_barrier + Vg_bias; % 게이트 
    
    Nd = [2e+20 ; 0 ; 2e+20]; % 도핑 농도 

    %% (.csv) 파일로부터 위치 정보를 정의 
    mesh_Generation();

    %% 상수 설정 - 글로벌 상수로 저장
    [ phi ] = initConstSet( Nd, Vgs );
                    
    %% Jacobian matrix 생성
    jbase = configueJbase();
    return;

    %% 초기 Poisson 방정식 정의
    [ phi, nn  ] = initPoisson2D( 100, jbase, phi);
    

% ### variable explorer code ###
%     phi = originVariable(1,1,phi);
%     nn = originVariable(1,1,nn);
%     global const_p
%     doping = const_p.doping;
%     doping = originVariable(1,1,doping);
    
%     global xmesh
%     x = xmesh.node;
%     global zmesh
%     z = zmesh.node';
    
sweep_mode = 1; % 1: Vds sweep, 2: Vgs sweep
if sweep_mode == 1
    %% bias sweep 설정 (Vds)
    bias = 0.5;     % 최종 bias 
    nVds = 51;      % bias point 개수
    Vds = linspace(0, bias, nVds);   % bias point 생성 
    Ids = Vds*0;    % 전류값을 저장할 변수 (Vds와 길이가 같아야 함)
    Vgs_delta = 0;  % Vgs를 sweep하지 않으므로 0
else
%     %% bias sweep 설정 (Vgs) - 이 기능을 위해서는 코드 수정이 필요 
%     bias = 0.5;     % 
%     nVgs = 51;      % 
%     Vgs_delta = linspace(0, bias, nVgs);     % simulation target 
%     Ids = Vgs_delta*0;  % 
% %     Vds = 0;        % Vds를 sweep하지 않으므로 0
end

%% self-consistent loop
for i = 1:51
    %% 해당 bias point에 대하여 미리 에너지 노드 설정과 각 노드에 대한 Fermi integral을 계산
    deltaE = 0.25e-3;    % 에너지 노드 간격 

    % valley 정보 : valley #1(l,t,t) #2(t,l,t) #3(t,t,l)
    
    % z방향 schrodinger 해석
    [Em_t, Vm_t] = mode_Confinement( 1, phi);   % valley #1, #2
    [Em_l, Vm_l] = mode_Confinement( 3, phi);   % valley #3

    % E_l 에너지(longitudinal energy) 범위 설정 
    E1 = min(min(Em_t(:,1)),min(Em_l(:,1)))-0.3;    % 시작점: Em 최소값보다 0.3 낮게 
    E2 = max(max(Em_t(:,2)),max(Em_l(:,5)))+0.3;    % 끝점: Em 최대값보다 0.3 높게 
    E1 = (round(E1/deltaE))*deltaE;     % 에너지 시작점을 deltaE 배수로 맞춤
    E2 = (round(E2/deltaE))*deltaE;     % 에너지 끝점을 deltaE 배수로 맞춤
    totalE = E1:deltaE:E2;  % 에너지 노드 생성

    nodeNum = 1000; % Fermi integral에 사용되는 E_y 에너지 적분 노드 갯수 
    % y방향 유효 질량 #2: m_y = m_l*m0, #1, #3: m_y = m_t*m0
    FF1_t = configue_Ffunction( 1, nodeNum, 0, totalE);  % #1, #3 / source
    FF2_t = configue_Ffunction( 1, nodeNum, Vds(i), totalE);  % #1, #3 / drain (bias)
    FF1_l = configue_Ffunction( 2, nodeNum, 0, totalE);     % #2 / source
    FF2_l = configue_Ffunction( 2, nodeNum, Vds(i), totalE);     % #2 / drain (bias)
    
    %% self-consistent loop 시작
    for j = 1:100
        phi_old = phi;
        tic;
       %% NEGF 해석
       % valley 정보 : #1(l,t,t) #2(t,l,t) #3(t,t,l)
       % z방향 schrodinger 해석
        [Em_t, Vm_t] = mode_Confinement( 1, phi);   % valley #1, #2 
        [Em_l, Vm_l] = mode_Confinement( 3, phi);   % valley #3 
        
        nn1 = phi*0; nn2 = nn1; nn3 = nn1; % nn1, nn2, nn3 저장변수 생성 
        for k = 1:2     % mode 갯수 2개 고려, 각 mode 에 대하여 계산되는 전자농도를 더함
            % find_Ffunction: Em영역에 따라 사용할 Fermi 함수 범위와 Energy 범위를 찾음
            % negf_Transport: NEGF를 해석하여 전자농도를 계산 (2곱하는 이유 - valley degeneracy)
            % valley #1(l,t,t)
            [E_v1, FF1_v1, FF2_v1] = find_Ffunction(Em_t, k, totalE, FF1_t, FF2_t); 
            nn1 = nn1 + 2*negf_Transport(1, Em_t, Vm_t, k, E_v1, FF1_v1, FF2_v1);
            % valley #2(t,l,t)
            [E_v2, FF1_v2, FF2_v2] = find_Ffunction(Em_t, k, totalE, FF1_l, FF2_l);
            nn2 = nn2 + 2*negf_Transport(2, Em_t, Vm_t, k, E_v2, FF1_v2, FF2_v2);            
        end
        for k = 1:5     % mode 갯수 5개 고려, 각 mode 에 대하여 계산되는 전자농도를 더함 
            % valley #3(t,t,l)
            [E_v3, FF1_v3, FF2_v3] = find_Ffunction(Em_l, k, totalE, FF1_t, FF2_t);
            nn3 = nn3 + 2*negf_Transport(3, Em_l, Vm_l, k, E_v3, FF1_v3, FF2_v3);
        end
        nn = nn1 + nn2 + nn3;   % 각 valley별 계산된 전자농도를 더함 
        
        %% Poisson 방정식 해석 
        % 최대 반복 횟수, jacobian 행렬, 포텐셜, 전자농도, Vgs변화값(기준은 초기 설정한 Vg_bias)
        [ phi, nn ] = nLinPoisson2D( 100, jbase, phi, nn, Vgs_delta);

        %% 수렴여부 체크하여 루프탈출여부 결정  
        stop(i,j) = max(max(abs(phi - phi_old)));
        disp(sprintf('[%d]self-consist loop[%d]-error: %d \n', i, j, stop(i,j)));
        if (stop(i,j) < 5e-4) || (j == 100)
            break;
        end
    end 
    %% 루프를 탈출하면 전류를 계산 (수렴하기 전 계산하는 것은 계산력 낭비)
    for k = 1:2     % mode 갯수 2개 고려, 각 mode 에 대하여 계산되는 전류값을 더함
        % valley #1(l,t,t)
        % negf_Current: NEGF를 해석하여 전류값을 계산 (2곱하는 이유 - valley degeneracy)
        [E_v1, FF1_v1, FF2_v1] = find_Ffunction(Em_t, k, totalE, FF1_t, FF2_t);
        Ids(i) = Ids(i) + 2*negf_Current(1, Em_t, k, E_v1, FF1_v1, FF2_v1);
        % valley #2(t,l,t)        
        [E_v2, FF1_v2, FF2_v2] = find_Ffunction(Em_t, k, totalE, FF1_l, FF2_l);
        Ids(i) = Ids(i) + 2*negf_Current(2, Em_t, k, E_v2, FF1_v2, FF2_);
    end
    for k = 1:5     % mode 갯수 5개 고려, 각 mode 에 대하여 계산되는 전류값을 더함
        % valley #3(t,t,l)
        [E_v3, FF1_v3, FF2_v3] = find_Ffunction(Em_l, k, totalE, FF1_t, FF2_t);
        Ids(i) = Ids(i) + 2*negf_Current(3, Em_l, k, E_v3, FF1_v3, FF2_v3);
    end
    
    %% 결과 저장 
    % 수렴된 결과를 각 bias point마다 저장. 
    if (Vds(i) - floor(Vds(i)*10)/10 < 1e-10) %% 0 V, 0.1 V, 0.2 V, 0.3 V, 0.4 V, 0.5 V
        % 0.1의 배수 point 마다 세부 계산 정보를 저장 
        % - 공간을 많이 차지하므로 원치 않으면 save함수만 남기고 주석처리하세요
        % 세부 계산 정보 : 각 valley, mode 별 spectral density A(x,E) (source/drain), 
        %                 transmission coefficient T(E), current Ids(E)
        % valley #1(l,t,t)
        [ v1_1_A1, v1_1_A2, v1_1_T, v1_1_Ids ] = negf_ShowVariables( 1, Em_t, 1, totalE, FF1_t, FF2_t );
        [ v1_2_A1, v1_2_A2, v1_2_T, v1_2_Ids ] = negf_ShowVariables( 1, Em_t, 2, totalE, FF1_t, FF2_t );
        % valley #2(t,l,t)
        [ v2_1_A1, v2_1_A2, v2_1_T, v2_1_Ids ] = negf_ShowVariables( 2, Em_t, 1, totalE, FF1_l, FF2_l );
        [ v2_2_A1, v2_2_A2, v2_2_T, v2_2_Ids ] = negf_ShowVariables( 2, Em_t, 2, totalE, FF1_l, FF2_l );
        % valley #3(t,t,l)
        [ v3_1_A1, v3_1_A2, v3_1_T, v3_1_Ids ] = negf_ShowVariables( 3, Em_l, 1, totalE, FF1_t, FF2_t );
        [ v3_2_A1, v3_2_A2, v3_2_T, v3_2_Ids ] = negf_ShowVariables( 3, Em_l, 2, totalE, FF1_t, FF2_t );
        [ v3_3_A1, v3_3_A2, v3_3_T, v3_3_Ids ] = negf_ShowVariables( 3, Em_l, 3, totalE, FF1_t, FF2_t );
        [ v3_4_A1, v3_4_A2, v3_4_T, v3_4_Ids ] = negf_ShowVariables( 3, Em_l, 4, totalE, FF1_t, FF2_t );
        [ v3_5_A1, v3_5_A2, v3_5_T, v3_5_Ids ] = negf_ShowVariables( 3, Em_l, 5, totalE, FF1_t, FF2_t );
        
        save(sprintf('result (Vg = 0.1V)\\result_%03d.mat',i)); % 저장
        % 저장 후 세부 정보 변수만 삭제 (메모리 점유율 낮추기 위함)
        clear v1_1_A1 v1_1_A2 v1_1_T v1_1_Ids 
        clear v1_2_A1 v1_2_A2 v1_2_T v1_2_Ids
        clear v2_1_A1 v2_1_A2 v2_1_T v2_1_Ids
        clear v2_2_A1 v2_2_A2 v2_2_T v2_2_Ids
        clear v3_1_A1 v3_1_A2 v3_1_T v3_1_Ids
        clear v3_2_A1 v3_2_A2 v3_2_T v3_2_Ids
        clear v3_3_A1 v3_3_A2 v3_3_T v3_3_Ids
        clear v3_4_A1 v3_4_A2 v3_4_T v3_4_Ids
        clear v3_5_A1 v3_5_A2 v3_5_T v3_5_Ids
    else % 수렴된 결과만 저장 
        save(sprintf('result (Vg = 0.1V)\\result_%03d.mat',i)); % 저장
    end
end 

return;





