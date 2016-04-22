function d = nearestApproachBrute(ra0, rb0, ra1, rb1, n)
% Compute the nearest approach of two line segments.
% Author: Jeff Comer <jcomer2@illinois.edu>

l0 = norm(rb0-ra0);
v0 = (rb0-ra0)/l0;
l1 = norm(rb1-ra1);
v1 = (rb1-ra1)/l1;
% Shift the origin to ra0.
r1 = ra1 - ra0;

dl0 = l0/n;
dl1 = l1/n;

d = inf;
for j=0:n
    for k=0:n
        d = min(d, norm(r1+v1*dl1*k - v0*dl0*j));
    end
end



