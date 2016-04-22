function r = moveDoubletTriplet(r)
% Make modified movement IV of Zhang et al. Biophys. J. 81 (2001) 1133-1443.
% Author: Jeff Comer <jcomer2@illinois.edu>
nNodes = length(r(:,1));

% Randomly choose doublets and triplets.
doubPick = ceil(rand*(nNodes-1));
tripPick = ceil(rand*(nNodes-2));
while doubPick >= tripPick-1 && doubPick <= tripPick+2
    doubPick = ceil(rand*(nNodes-1));
    tripPick = ceil(rand*(nNodes-2));
end

d0 = r(doubPick,:);
d10 = r(doubPick+1,:) - d0;
t0 = r(tripPick,:);
t20 = r(tripPick+2,:) - t0;
t10 = r(tripPick+1,:) - t0;

% z is the progress that t1 makes from t0 to t2.
z = dot(t10,t20)/dot(t20,t20);
s = sqrt(dot(t10,t10) - z^2);

% Choose a new basis along d10 with a randomly chosen azimuthal angle.
phi = 2.0*pi*rand;
a = [0 0 1];
ez = d10/norm(d10);
ex = a - dot(a,ez)*ez;
ex = ex/norm(ex);
es = ex*cos(phi) + cross(ez,ex)*sin(phi);

% Compute the position of the new node.
pos = d0 + z*ez + s*es;

% Reform the position list.
if doubPick < tripPick
    r = [r(1:doubPick,:); pos; r(doubPick+1:tripPick,:); r(tripPick+2:end,:)];
else
    r = [r(1:tripPick,:); r(tripPick+2:doubPick,:); pos; r(doubPick+1:end,:)];
end



