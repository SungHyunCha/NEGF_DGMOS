  function [phi_new, nn_new]  = nLinPoisson2D(iterNum, jbase, nxz, ni, phi_new, boundary, doping, ...
                    x, x_dlt, x_int, z, z_dlt, z_int, eps0, Vt, q)
% function [phi] = nLinPoisson2D(iterNum, jbase, nxz, ni, phi, boundary, Vgs, doping, ...
%                     x_dlt, x_int, z_dlt, z_int, eps0, Vt, q)

% set index x-direction
nx     = size(phi_new,1);
x_int1 = x_int{1};
x_int2 = x_int{2};

% set index z-direction
nz     = size(phi_new,2);
z_int1 = z_int{1};
z_int2 = z_int{2};

BC_index = find(boundary ~= 0);

% ### Vgs sweep only ### 
% boundary(BC_index) = boundary(BC_index) + Vgs;

oldphi = phi_new;

% doping(:,z_int1) = 0;
% doping(:,z_int2) = 0;

% for test
for i = 1:iterNum
    %% #1. residual vector 
    % calculating charge density 
    nn_new = nxz.*exp(+(phi_new-oldphi)/Vt); 
    pp_new = ni.*exp(-phi_new/Vt);
    g  = q*(pp_new - nn_new + doping);  
    g_log=g;
    
    % interface exception (:g)
    g(x_int1,:) = bsxfun( @times, g(x_int1,:), x_dlt(x_int1-1)./(x_dlt(x_int1-1)+x_dlt(x_int1)) ) ...
                + bsxfun( @times, g(x_int2,:), x_dlt(x_int1)  ./(x_dlt(x_int1-1)+x_dlt(x_int1)) );
    g(:,z_int1) = bsxfun( @times, g(:,z_int1), (z_dlt(z_int1-1)./(z_dlt(z_int1-1)+z_dlt(z_int1)))' ) ...
                + bsxfun( @times, g(:,z_int2), (z_dlt(z_int1)  ./(z_dlt(z_int1-1)+z_dlt(z_int1)))' );
    g(x_int2,:) = 0;    
    g(:,z_int2) = 0;
    
    % boundary exception 
    g(BC_index) = 0;
    
    r = ( jbase*matrixToVector(phi_new, nx, nz) - matrixToVector(boundary, nx, nz) )...
        *eps0*1e+9^2;
    R = r + matrixToVector(g, nx, nz);

    %% #2. jacobian matrix 
    h = -(q/Vt)*(pp_new + nn_new);
    % interface exception (:h)
    h(x_int1,:) = bsxfun( @times, h(x_int1,:), x_dlt(x_int1-1)./(x_dlt(x_int1-1)+x_dlt(x_int1)) ) ...
                + bsxfun( @times, h(x_int2,:), x_dlt(x_int1)  ./(x_dlt(x_int1-1)+x_dlt(x_int1)) );
    h(:,z_int1) = bsxfun( @times, h(:,z_int1), (z_dlt(z_int1-1)./(z_dlt(z_int1-1)+z_dlt(z_int1)))' ) ...
                + bsxfun( @times, h(:,z_int2), (z_dlt(z_int1)  ./(z_dlt(z_int1-1)+z_dlt(z_int1)))' );
    h(x_int2,:) = 0;
    h(:,z_int2) = 0;
    
    % boundary exception 
    h(BC_index) = 0;
    
    j = jbase*(eps0)*1e+9^2;
    J = j + diag( matrixToVector(h, nx, nz) );
    
    %% #3. finding negative delta vector 
    invJ = inv(J);
    dphi = -invJ*R;
    dphiMat = vectorToMatrix(dphi,nx,nz);
    
    phi_new = phi_new + dphiMat;
    dphiVec = matrixToVector(dphiMat, nx, nz);
    stop = full(max(abs(dphiVec)));
%     stop_log1(i,1) = stop;

    disp(sprintf('nonLinear Poisson trial[%d]-error: %d \n', i, stop))
    if stop < 1e-12
%         save nLinPoisson
        break;
    end
end

end 
% 
