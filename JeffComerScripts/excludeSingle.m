function [reject collideA collideB] = excludeSingle(r, diam)
% Compute the nearest approach of two line segments.
% Author: Jeff Comer <jcomer2@illinois.edu>
nNodes = length(r(:,1));
collideA = [0 0 0];
collideB = [0 0 0];

reject = 0;
for j=1:nNodes-2
    for k=j+2:nNodes-1
        d = nearestApproachBrute(r(j,:),r(j+1,:),r(k,:),r(k+1,:),20);
        
        if d < diam
            reject = 1;
            collideA = r(j:j+1,:);
            collideB = r(k:k+1,:);
        end
    end
    if reject, break, end
end



