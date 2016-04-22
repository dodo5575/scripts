clear all;

%inDir = 'tilt';
%run='tilt';
%nameList = {'tilt10cM_at_neg' 'tilt10cM_at_pos' 'tilt10cM_gc_neg' 'tilt10cM_gc_pos'};

inDir = 'pmf_basepair';
run='best140cM_all';
nameList = {'best140cM_at_neg' 'best140cM_at_pos' 'best140cM_gc_neg' 'best140cM_gc_pos'...
    'best140cM1_at_neg' 'best140cM1_at_pos' 'best140cM1_gc_neg' 'best140cM1_gc_pos'...
'best140cM*_at_neg' 'best140cM*_at_pos' 'best140cM*_gc_neg' 'best140cM*_gc_pos'};

timestep = 2e-5; % in ns
outPeriod = 2000;
cutIndex = 10;
dt = outPeriod*timestep;

fprintf('\n%s\n', run);
out = fopen(sprintf('current_%s.dat', run), 'w');
for j=1:length(nameList)
    data = [];
    glob = sprintf('%s/*%s*.curr', inDir, nameList{j});
    fileList = dir(glob);
    nFiles = length(fileList);
    %fprintf('Found %d files.\n', nFiles);

    if nFiles == 0
        continue
    end
        
    for f=1:nFiles
        fileName = sprintf('%s/%s', inDir, fileList(f).name);
        %fprintf('file %s\n', fileName);
        if exist(fileName, 'file')
            data0 = dlmread(fileName, ' ');
            data = [data; data0(cutIndex:end,:)];
        else
            fprintf('%s does not exist\n', fileName);
        end
    end
    
    n = length(data(:,1));
    meanCurr = 1e3*mean(data(:,2));
    stdCurr = 1e3*std(data(:,2));
    errCurr = stdCurr/sqrt(n);
    
    fprintf(out, '%s  \t\t%+.2f\t\t%.2f\t\t%.1f\n', nameList{j}, meanCurr, errCurr, n*dt);
    disp(sprintf('%s  \t\t%+.2f\t\t%.2f\t\t%.1f', nameList{j}, meanCurr, errCurr, n*dt));
end

fclose(out);
