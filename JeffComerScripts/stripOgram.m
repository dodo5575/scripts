clear all;

inSet = [0:24 300:315];
inName = 'middling';
inPrefix = sprintf('%s/strip_%s_pmf', inName, inName);
inSuffix = '.dat';

depth = zeros(length(inSet), 1);

for k=1:length(inSet)
    % Read the profile.    
    profFile = sprintf('%s%d%s', inPrefix, inSet(k), inSuffix);
    data = dlmread(profFile, ' ');
    profZ = data(:,1);
    profU = data(:,2);
    
    depth(k) = min(profU);
end


% Histogram the data.
[n f] = hist(depth);
plot(f, n, 'b-');

outName = sprintf('hist_%s.dat', inName);
dlmwrite(outName, [f' n'], ' ');