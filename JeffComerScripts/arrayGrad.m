function [ax ay] = arrayGrad(a, nx, h)
% Compute the gradient of serialized double array (x fast, y slow).
% A second order finite difference and periodic boundary conditions are
% used.
% Author: Jeff Comer <jcomer2@illinois.edu>
n = length(a);
ny = ceil(n/nx);
ax = zeros(size(a));
ay = zeros(size(a));

% Do the internal nodes.
for iy=1:ny-2
    for ix=1:nx-2
        j = 1 + ix + nx*iy;
        
        ax(j) = 0.5*(a(j+1)-a(j-1))/h;
        ay(j) = 0.5*(a(j+nx)-a(j-nx))/h;
    end
end

% Do the top and bottom.
for ix=1:nx-2
    jt = 1 + ix;
    jb = 1 + ix + nx*(ny-1);

    ax(jt) = 0.5*(a(jt+1)-a(jt-1))/h;
    ax(jb) = 0.5*(a(jb+1)-a(jb-1))/h;
    
    ay(jt) = 0.5*(a(jt+nx)-a(jb))/h;
    ay(jb) = 0.5*(a(jt)-a(jb-nx))/h;
end

% Do the left and right.
for iy=1:ny-2
    jl = 1 + nx*iy;
    jr = 1 + nx-1 + nx*iy;
    
    ax(jl) = 0.5*(a(jl+1)-a(jr))/h;
    ax(jr) = 0.5*(a(jl)-a(jr-1))/h;
    
    ay(jl) = 0.5*(a(jl+nx)-a(jl-nx))/h;
    ay(jr) = 0.5*(a(jr+nx)-a(jr-nx))/h;
end

% Do the corners.
j00 = 1;
j10 = 1 + nx;
j01 = 1 + nx*(ny-1);
j11 = nx*ny;

ax(j00) = 0.5*(a(j00+1)-a(j10))/h;
ax(j10) = 0.5*(a(j00)-a(j10-1))/h;
ax(j01) = 0.5*(a(j01+1)-a(j11))/h;
ax(j11) = 0.5*(a(j01)-a(j11-1))/h;

ay(j00) = 0.5*(a(j00+nx)-a(j01))/h;
ay(j10) = 0.5*(a(j10+nx)-a(j11))/h;
ay(j01) = 0.5*(a(j00)-a(j01-nx))/h;
ay(j11) = 0.5*(a(j10)-a(j11-nx))/h;



