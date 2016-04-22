function m = arrayDeserialize(a, nx)
% Converts serialized double array (x fast, y slow) into a matrix.
% Author: Jeff Comer <jcomer2@illinois.edu>
n = length(a);
ny = ceil(n/nx);

if round(ny) ~= ny || ny < 1
    error('arrayDeserialize:invalidDim', 'Array dimension is invalid!');
end
m = zeros(nx,ny);

for iy=1:ny
    j0 = 1 + (iy-1)*nx;
    j1 = iy*nx;
    m(:,iy) = a(j0:j1);
end





