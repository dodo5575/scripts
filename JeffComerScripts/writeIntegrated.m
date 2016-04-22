clear all;
% Author: Jeff Comer <jcomer2@illinois.edu>
glob = '*.vcorr';
corrFactor = 1e6;

fileList = dir(glob);
nFiles = length(fileList);
fprintf('Found %d files.\n', nFiles);

for j=1:nFiles
  fileName = sprintf('%s', fileList(j).name);
  
  data = dlmread(fileName, ' ');
  [t diff] = integrateDiffusivity(data(:,1), data(:,2));
  diff = diff*corrFactor;
  
  outFile = sprintf('%s.int', fileName);
  dlmwrite(outFile, [t diff], ' ');
  fprintf('Wrote "%s".\n', outFile);
end
