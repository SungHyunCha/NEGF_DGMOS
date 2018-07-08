function [ phi ] = initPotential( phi, doping, ni, Vt )
%% 함수설명 : 도핑농도에 따른 초기 포텐셜 값을 입력합니다. 
% 각 파라미터는 다음과 같습니다. 
% phi : 미리 정의된 포텐셜 저장 공간 
% doping : 위치에 따른 도핑 농도 
% ni : intrinsic 전자 농도 
% Vt : Thermal voltage 

%% 포텐셜 입력 
% intrinsic에 대하여는 0 으로 포텐셜 입력 

pp_index = find(doping < 0);    % p-type 에 대하여 
phi(pp_index) = -(Vt)*log(-doping(pp_index)./ni(pp_index)); % - 포텐셜 입력 

nn_index = find(doping > 0);    % n-type 에 대하여 
phi(nn_index) = +(Vt)*log(+doping(nn_index)./ni(nn_index)); % + 포텐셜 입력 
end
