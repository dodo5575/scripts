function force = interpolateForce(data, grid, delta, origin, pos)
% Cubic interpolation of the derivative of a grid
fx = interpGrad(data, grid, delta, origin, pos, 1);
fy = interpGrad(data, grid, delta, origin, pos, 2);
fz = interpGrad(data, grid, delta, origin, pos, 3);
f = [fx; fy; fz];
force = inv(delta')*f;

% The derivative is taken along direction dir.
function f = interpGrad(data, grid, delta, origin, pos, direc)
% Find the home node.
l = inv(delta)*(pos-origin);
home = floor(l);
%f = 0.0;
%if any(home < 1) || any(home >= grid-2) return, end

% Choose the direction.
dirY = mod(direc, 3) + 1;
dirZ = mod(direc + 1, 3) + 1;
dir = [direc dirY dirZ];
jump = [grid(3)*grid(2) grid(3) 1];
jump = jump(dir);
me = home(dir);
g = grid(dir);

% Find the interpolation coordinates.
w = l - home;
w = w(dir);

% Find the values at the neighbors.
g1 = zeros(4,4,4);
for ix=1:4
  for iy=1:4
    for iz=1:4
      % Wrap around the periodic boundaries.
      jx = mod(ix-2 + me(1), g(1));  
      jy = mod(iy-2 + me(2), g(2));  
      jz = mod(iz-2 + me(3), g(3));  
        
      ind = jz*jump(3) + jy*jump(2) + jx*jump(1);
      g1(ix,iy,iz) = data(ind+1);
    end
  end
end

% Mix along x, taking the derivative.
g2 = zeros(4,4);
for iy=1:4
  for iz=1:4
    a3 = 0.5*(-g1(1,iy,iz) + 3*g1(2,iy,iz) - 3*g1(3,iy,iz) + g1(4,iy,iz));
    a2 = 0.5*(2*g1(1,iy,iz) - 5*g1(2,iy,iz) + 4*g1(3,iy,iz) - g1(4,iy,iz));
    a1 = 0.5*(-g1(1,iy,iz) + g1(3,iy,iz));
    a0 = g1(2,iy,iz);

    g2(iy,iz) = 3.0*a3*w(1)*w(1) + 2.0*a2*w(1) + a1;
    %g2(iy,iz) = a3*w(1)*w(1)*w(1) + a2*w(1)*w(1) + a1*w(1) + a0;
  end
end

% Mix along y.
g3 = zeros(4,1);
for iz=1:4
   a3 = 0.5*(-g2(1,iz) + 3*g2(2,iz) - 3*g2(3,iz) + g2(4,iz));
   a2 = 0.5*(2*g2(1,iz) - 5*g2(2,iz) + 4*g2(3,iz) - g2(4,iz));
   a1 = 0.5*(-g2(1,iz) + g2(3,iz));
   a0 = g2(2,iz);
   
   g3(iz) = a3*w(2)*w(2)*w(2) + a2*w(2)*w(2) + a1*w(2) + a0;
end

% Mix along z.
a3 = 0.5*(-g3(1) + 3*g3(2) - 3*g3(3) + g3(4));
a2 = 0.5*(2*g3(1) - 5*g3(2) + 4*g3(3) - g3(4));
a1 = 0.5*(-g3(1) + g3(3));
a0 = g3(2);

f = -(a3*w(3)*w(3)*w(3) + a2*w(3)*w(3) + a1*w(3) + a0);
