function [ E, FF1, FF2 ] = find_Ffunction( Em, k_count, totalE, FF1, FF2 )
%% 함수설명 : total 에너지에서 정의된 FF1, FF2 중 계산에 필요한 영역을 찾아서 반환
% 각 파라미터는 다음과 같습니다. 
% Em : subband minimum
% k_count : 해석할 mode 번호 
% totalE : 시뮬레이션에서 고려되는 모든 영역의 에너지 
% FF1 : Source단에서의 Fermi 함수값, E index에 대응
% FF2 : Drain단에서의 Fermi 함수값, E index에 대응

%% 계산에 필요한 영역을 찾아서 반환 
deltaE = totalE(2) - totalE(1); % 에너지 노드 간격 
Em = Em(:,k_count); % 주어진 mode 번호에 대한 subband minimum 선택
E1 = min(Em)-0.15;  % subband minimum보다 0.15 eV 낮은 에너지부터 
E2 = max(Em)+0.25;  % subband maximum보다 0.25 eV 높은 에너지까지만 고려 
E1 = (round(E1/deltaE))*deltaE;     % 영역 첫 지점을 에너지 노드 간격의 배수로 맞춤 
E2 = (round(E2/deltaE))*deltaE;     % 영역 끝 지점을 에너지 노드 간격의 배수로 맞춤 
ind1 = find(abs(totalE - E1) < deltaE/5, 1, 'first');   % 영역 첫 지점의 인덱스를 총 에너지 영역에서 가져옴
ind2 = find(abs(totalE - E2) < deltaE/5, 1, 'first');   % 영역 끝 지점의 인덱스를 총 에너지 영역에서 가져옴

% 빈 영역의 인덱스를 참조하였을 경우에 대한 예외처리 
i = isempty(ind1);
j = isempty(ind2);
if (i == 1) || (j == 1) % 둘 중에 하나라도 비어있으면 처리 시작 
    disp('failure to find sub-energy array in total-energy region');    % 계산에 문제가 있음을 알림 
    if (i == 1) % 영역 첫 지점을 못 찾았으면 
        ind1 = 1;   % 첫 번째 인덱스로 강제지정 
    end
    if (j == 1) % 끝 지점을 못찾았으면 
        ind2 = size(totalE,2);  % 마지막 인덱스로 강제지정 
    end
end

% 첫 지점과 끝 지점까지의 인덱스를 이용, 필요한 부분을 찾음 
% Fermi함수의 길이가 하나 짧은 이유는 노드 간 중점에서의 결과값이기 때문. 
% NEGF 처리 시 에너지 노드에 대하여 같은 처리를 수행하므로 문제없음. 
E = totalE(ind1:ind2);      % 에너지 영역 
FF1 = FF1(ind1:ind2-1);     % Fermi 함수(source)
FF2 = FF2(ind1:ind2-1);     % Fermi 함수(drain)
end