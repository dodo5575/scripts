function b = arrayBlur3d(a, nx, ny, nz)
% Compute the gradient of serialized triple array (x fast, y medium, z slow).
% A second order finite difference and periodic boundary conditions are
% used.
% jcomer2@uiuc.edu
b = zeros(size(a));

for iz=0:nz-1
    for iy=0:ny-1
        for ix=0:nx-1
            k = getNeighbors([ix iy iz], [nx ny nz]);
            
            b(k(1)) = mean(a(k(2:end)));
        end
    end
end
