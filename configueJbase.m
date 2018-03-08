function [ jbase ] = configueJbase(eps_r, boundary, x_idx, x_int, x_dlt, z_idx, z_int, z_dlt)
%% index setting
% set index x-direction
nx = x_idx{end};
x_int1 = x_int{1};
x_int2 = x_int{2};

% set index z-direction
nz = z_idx{end};
z_int1 = z_int{1};
z_int2 = z_int{2};

% configue delta profile 
x_dlt(x_int1) = x_dlt(x_int2);

% configue delta profile 
z_dlt(z_int1) = z_dlt(z_int2);

% set metal line 
boundary = boundary(:,:,1);

%% jacobi matrix setting 
% jacobi matrix definition 
jbase = zeros(nx*nz, nx*nz, 2);

% config left eps_r profile 
eps_x = eps_r;
eps_x(:,z_int1) = bsxfun(@times, eps_r(:,z_int1), (z_dlt(z_int1-1)./(z_dlt(z_int1-1)+z_dlt(z_int2)))' ) ...
                + bsxfun(@times, eps_r(:,z_int2), (z_dlt(z_int2)  ./(z_dlt(z_int1-1)+z_dlt(z_int2)))' );
            
% eps_x(:,z_int1(1)) = eps_r(:,z_int1(1));    % oxide 
% eps_x(:,z_int1(2)) = eps_r(:,z_int2(2));    % oxide 

% eps_x(:,z_int1(1)) = eps_r(:,z_int2(1));    % silicon
% eps_x(:,z_int1(2)) = eps_r(:,z_int1(2));    % silicon 

            
eps_left        = eps_x(2:end,:);
eps_left(end,:) = eps_x(end,:)*2;

% configue right epsilon profile
eps_right           = eps_x(1:end-1,:);
eps_right(x_int1,:) = eps_x(x_int2,:);  % x_int2 govern x-interface at right case
eps_right(1,:)      = eps_x(1,:)*2;

% configue lower epsilon profile 
eps_z = eps_r;
eps_z(x_int1,:) = bsxfun(@times, eps_r(x_int1,:), x_dlt(x_int1-1)./(x_dlt(x_int1-1)+x_dlt(x_int2)) ) ...
                + bsxfun(@times, eps_r(x_int2,:), x_dlt(x_int2)  ./(x_dlt(x_int1-1)+x_dlt(x_int2)) );
            
eps_lower        = eps_z(:,2:end);
eps_lower(:,end) = eps_z(:,end)*2;

% configue upper epsilon profile 
eps_upper           = eps_z(:,1:end-1);
eps_upper(:,z_int1) = eps_z(:,z_int2);  % z_int2 govern z-interface at upper case
eps_upper(:,1)      = eps_z(:,1)*2;

%% base matrix configue 
% jacobi left basis definition 
delta_n = x_dlt; 
delta_p = [x_dlt(2:end) ; x_dlt(end)];
j_left  = zeros(nx, nz);
j_left(2:end,:)  = bsxfun(@rdivide, eps_left, delta_n.*(delta_p+delta_n)/2 );
% j_left(2:end,:)  = bsxfun(@rdivide, eps_left, 1);
j_cent1 = -j_left;

% jacobi right basis definition 
delta_n = [x_dlt(1) ; x_dlt(1:end-1)]; 
delta_p = x_dlt;
j_right  = zeros(nx, nz);
j_right(1:end-1,:)  = bsxfun(@rdivide, eps_right, delta_p.*(delta_p+delta_n)/2 );
% j_right(1:end-1,:)  = bsxfun(@rdivide, eps_right, 1);
j_cent2 = -j_right;

% jacobi lower basis definition 
delta_n = z_dlt;
delta_p = [z_dlt(2:end) ; z_dlt(end)];
j_lower  = zeros(nx, nz);
j_lower(:,2:end) = bsxfun(@rdivide, eps_lower, (delta_n.*(delta_p+delta_n)/2)' );
% j_lower(:,2:end) = bsxfun(@rdivide, eps_lower, 1);
j_cent3 = -j_lower;

% jacobi upper basis definition 
delta_n = [z_dlt(1) ; z_dlt(1:end-1)];
delta_p = z_dlt;
j_upper  = zeros(nx, nz);
j_upper(:,1:end-1) = bsxfun(@rdivide, eps_upper, (delta_p.*(delta_p+delta_n)/2)' );
% j_upper(:,1:end-1) = bsxfun(@rdivide, eps_upper, 1);
j_cent4 = -j_upper;

% total jacobi configuration (x-axis)
j_cent1 = matrixToVector(j_cent1,nx,nz);  % center1
j_left  = matrixToVector(j_left, nx,nz);  % left 
j_right = matrixToVector(j_right,nx,nz);  % right 
j_cent2 = matrixToVector(j_cent2,nx,nz);  % center2 

jbase(2:end, 2:end,     1) = jbase(2:end, 2:end,     1) + diag(j_cent1(2:end));
jbase(2:end, 1:end-1,   1) = jbase(2:end, 1:end-1,   1) + diag(j_left (2:end));
jbase(1:end-1, 2:end,   1) = jbase(1:end-1, 2:end,   1) + diag(j_right(1:end-1));
jbase(1:end-1, 1:end-1, 1) = jbase(1:end-1, 1:end-1, 1) + diag(j_cent2(1:end-1));

% total jacobi configuration (z-axis)
j_cent3 = matrixToVector(j_cent3,nx,nz);  % center3
j_lower = matrixToVector(j_lower,nx,nz);  % lower
j_upper = matrixToVector(j_upper,nx,nz);  % upper 
j_cent4 = matrixToVector(j_cent4,nx,nz);  % center4

jbase(1+nx:end, 1+nx:end, 2) = jbase(1+nx:end, 1+nx:end, 2) + diag(j_cent3(1+nx:end));
jbase(1+nx:end, 1:end-nx, 2) = jbase(1+nx:end, 1:end-nx, 2) + diag(j_lower(1+nx:end));
jbase(1:end-nx, 1+nx:end, 2) = jbase(1:end-nx, 1+nx:end, 2) + diag(j_upper(1:end-nx));
jbase(1:end-nx, 1:end-nx, 2) = jbase(1:end-nx, 1:end-nx, 2) + diag(j_cent4(1:end-nx));

jbase = jbase(:,:,1) + jbase(:,:,2);


%% interface process 
% interface process (x-axis)
for i = 1:size(x_int1,1)
    for j = 0:nz-1
        x1 = x_int1(i) + j*nx;
        x2 = x_int2(i) + j*nx;
        jbase(x1,x2+1) = jbase(x1,x2);
        jbase(x1,x2  ) = 0;

        jbase(x2,:)  = 0;
        jbase(x2,x2) = 1;
        jbase(x2,x1) = -1;
    end 
end

% interface procees (z-axis)
for i = 1:nx
    for j = 1:size(z_int1,1)
        z1 = i+(z_int1(j)-1)*nx;
        z2 = i+(z_int2(j)-1)*nx;
        jbase(z1,z2+nx) = jbase(z1,z2);
        jbase(z1,z2   ) = 0;

        jbase(z2,:)  = 0;
        jbase(z2,z2) = 1;
        jbase(z2,z1) = -1;
    end 
end

% interface process (xz-intersection)
for i = 1:size(x_int1,1)
   for j  = 1:size(z_int1,1)
       xz1 = x_int1(i) + (z_int1(j)-1)*nx;
       xz2 = x_int2(i) + (z_int2(j)-1)*nx;
       
       jbase(xz2, :) = 0;
       jbase(xz2, xz2) = 1; 
       jbase(xz2, xz1) = -1;
   end
end

max(sum(jbase,2))

% boundary setting (Dirichlet boundary) 
jbase_BC = eye(nx*nz);
BC_index = find(boundary ~= 0);
jbase(BC_index,:) = jbase_BC(BC_index,:);    % bottom gate


end
























