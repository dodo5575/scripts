function r = moveStretch(r,distMax)
% Make modified movement III of Zhang et al. Biophys. J. 81 (2001) 1133-1443.
% Author: Jeff Comer <jcomer2@illinois.edu>
nNodes = length(r(:,1));

% Randomly choose a segment.
pick = ceil(rand*(nNodes-1));
d = r(pick+1,:)-r(pick,:);
% Make a random displacement along the direction between the segments.
dr = d/norm(d)*distMax*(2*rand-1);

for j=pick+1:nNodes
    r(j,:) = r(j,:) + dr;
end



