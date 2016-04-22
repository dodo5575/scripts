function r = moveRotation(r,angleMax)
% Make modified movement II of Zhang et al. Biophys. J. 81 (2001) 1133-1443.
% Author: Jeff Comer <jcomer2@illinois.edu>
nNodes = length(r(:,1));

% Randomly choose a rotation center, magnitude, and axis.
pick = floor(rand*nNodes);
angle = angleMax*rand;
theta = 2*pi*rand;
phi = acos(2*rand-1);
a = [sin(phi)*cos(theta) sin(phi)*sin(theta) cos(phi)];

if pick == 0
    r0 = [0 0 0];
else
    r0 = r(pick,:);
end

for j=pick+1:nNodes
    rj = r(j,:)-r0;
    z = a*dot(a,rj);
    r(j,:) = (rj-z)*cos(angle) + cross(a,rj)*sin(angle) + z + r0;
end
