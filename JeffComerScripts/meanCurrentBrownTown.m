clear all;

inName = '';
dir = 'high';
run = 'high';
outPeriod = 500;
timestep = 2e-5; % in ns
rawPrefix = 'raw/rawcurr';
%nameList = {'open1_vol' 'open1_vol_spots' 'open1_spots' 'open1_third'
%'dna_vol' 'dna_vol_spots' 'dna_spots' 'dna_third'};
nameList = {'open_spots' 'open_nospots' 'dna_spots' 'dna_nospots'};
ionList = {''};
sims = 20;
cutIndex = 100;

simList = cell(sims,1);
for j=1:sims
    simList{j} = num2str(j-1);
end
for j=1:length(nameList)
    currList{j} = simList;
end

outName = sprintf('mean_%s_browntown.txt', run);
out = fopen(outName, 'w');

fprintf('\n%s\n', run);
fprintf(out, '%s\n', run);
fprintf(out, 'system\t\tcurrent (pA)\t\tse (pA)\t\ttime (ns)');
fprintf('system\t\tcurrent (pA)\t\tse (pA)\t\ttime (ns)');
dcdFreq = outPeriod*timestep;

for k=1:length(ionList)
  fprintf(out, '\nion %s\n', ionList{k});
  disp(sprintf('\nion %s', ionList{k}));
  for j=1:length(nameList)
    data = [];
    
    currSet = currList{j};
    %currSet = currSet(2:end);
    
    for c=1:length(currSet)
      fileName = sprintf('%s/%s.brown.%s.curr', dir, nameList{j}, currSet{c});
      if exist(fileName, 'file')
        data0 = dlmread(fileName, ' ');
        data = [data; data0(cutIndex:end,:)];
      end
    end
    
    n = length(data(:,1));
    meanCurr = 1e3*mean(data(:,2));
    stdCurr = 1e3*std(data(:,2));
    errCurr = stdCurr/sqrt(n);
    
    rawFile = sprintf('%s%s_%s_%s.dat', rawPrefix, ionList{k}, nameList{j}, dir);
    fprintf(out, '%s  \t\t%+.2f\t\t%.2f\t\t%.1f\n', nameList{j}, meanCurr, errCurr, n*dcdFreq);
    disp(sprintf('%s  \t\t%+.2f\t\t%.2f\t\t%.1f', nameList{j}, meanCurr, errCurr, n*dcdFreq));
  end
end

fclose(out);
