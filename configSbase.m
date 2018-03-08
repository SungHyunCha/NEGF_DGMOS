function [ sbase1, sbase2 ] = configSbase( nz )
%% configue Schrodinger base matrix 
sbase_first_row = zeros(nz,1);
sbase_first_col = zeros(1,nz);
sbase_first_row(1) = -1; sbase_first_row(2) = 1;
sbase_first_col(1) = -1; sbase_first_col(2) = 0; 
sbase1 = sparse(toeplitz(sbase_first_row, sbase_first_col));    % left 
sbase2 = sbase1';   % right 
end

