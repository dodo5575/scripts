function r = randomConformationRigid(n, len)
% Author: Jeff Comer <jcomer2@illinois.edu>

r = zeros(n,3);
for j=2:n
    % Pick angles uniform on a sphere.
    theta = 2*pi*rand;
    phi = acos(2*rand-1);
    dr = len*[sin(phi)*cos(theta) sin(phi)*sin(theta) cos(phi)];
    
    r(j,:) = r(j-1,:) + dr;
end



