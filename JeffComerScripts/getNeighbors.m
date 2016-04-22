function neigh = getNeighbors(p, n)
% Given a position p = [ix iy iz] and n = [nx ny nz], we find the nearest
% neighbors of the home node with x fast, y medium, z slow.
% The first neighbor is the home node.
% Note we take p to be zero-based.
% jcomer2@uiuc.edu

% Get the eight neighbors with wrapping.
n0 = zeros(1,3);
n1 = zeros(1,3);
for k=1:3
    if p(k) == 0
        n0(k) = n(k)-1;
    else
        n0(k) = p(k)-1;
    end
    
    if p(k) == n(k)-1
        n1(k) = 0;
    else
        n1(k) = p(k)+1;
    end
end

% Make a list of indices.
nx = n(1);
nxny = n(1)*n(2);
j = 1 + p(1) + p(2)*nx + p(3)*nxny;
jx0 = 1 + n0(1) + p(2)*nx + p(3)*nxny;
jx1 = 1 + n1(1) + p(2)*nx + p(3)*nxny;
jy0 = 1 + p(1) + n0(2)*nx + p(3)*nxny;
jy1 = 1 + p(1) + n1(2)*nx + p(3)*nxny;
jz0 = 1 + p(1) + p(2)*nx + n0(3)*nxny;
jz1 = 1 + p(1) + p(2)*nx + n1(3)*nxny;
neigh = [j jx0 jx1 jy0 jy1 jz0 jz1];
