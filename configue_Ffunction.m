function [ F ] = configue_Ffunction(valleyNum, nodeNum, Vds, E)
%% 함수설명 : homogeneous 방향으로 적분된 Fermi 함수를 계산합니다.
% 각 파라미터는 다음과 같습니다. 
% valleyNum : y방향으로의 valley 번호 (#2: m_l, #1 & #3 : m_t)
% nodeNum : 적분을 위해 homogeneous 방향 에너지에 할당할 노드 개수 (높을수록 정확)
% Vds : bias에 의해 발생하는 Fermi 에너지의 offset [V]
% E : 해석에 사용할 longitudinal 에너지 [eV]

%% 글로벌 상수를 불러옵니다. 
if (3 < valleyNum ) || (valleyNum < 1)
    disp(sprintf('option : out of range!'));
    return;
end
global mass;    % valley에 따른 전자 유효 질량 
m_y = mass.m_y(valleyNum);

% NEGF에서 사용할 구분구적법을 위해 에너지 노드를 중간값으로 변환 
dE = E(2) - E(1);   
E_l = E(1:end-1) + dE/2; % 변환된 longitudinal 에너지

% 특정 값을 가지는 longitudinal 에너지에 대하여 homogeneous 방향으로 적분된 Fermi 함수를 계산
F = Fhalfinv( -Vds -E_l, nodeNum, m_y);

end

