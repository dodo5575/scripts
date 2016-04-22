% Make a color map from a list colors.
% Author: Jeff Comer <jcomer@illinois.edu>
function map = makeColorMap(colorList, samples)

n = length(colorList(:,1))-1;
delta = samples/n;

map = zeros(samples, 3);
for j=1:samples
    home = floor((j-1)/delta);
    next = home + 1;
    homePos = home*delta;
    alpha = (j-1-homePos)/delta;

    % Mix the colors.
    map(j,:) = colorList(home+1,:)*(1-alpha) + colorList(next+1,:)*alpha;
    %map(j,:) = [home next alpha];
end
