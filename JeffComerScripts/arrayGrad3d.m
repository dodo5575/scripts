function [ax ay az] = arrayGrad3d(a, nx, ny, nz, h)
% Compute the gradient of serialized triple array (x fast, y medium, z slow).
% A second order finite difference and periodic boundary conditions are
% used.
% jcomer2@uiuc.edu
ax = zeros(size(a));
ay = zeros(size(a));
az = zeros(size(a));

for iz=0:nz-1
    for iy=0:ny-1
        for ix=0:nx-1
            k = getNeighbors([ix iy iz], [nx ny nz]);
            
            ax(k(1)) = 0.5*(a(k(3))-a(k(2)))/h;
            ay(k(1)) = 0.5*(a(k(5))-a(k(4)))/h;
            az(k(1)) = 0.5*(a(k(7))-a(k(6)))/h;
        end
    end
end
