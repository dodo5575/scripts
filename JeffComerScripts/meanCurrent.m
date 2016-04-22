clear all;

%inName = '';
%dir = 'current_all';
inName = 'thru28_';
dir = 'current_through';
outName = sprintf('mean_%s.txt', dir);
rawPrefix = 'raw/rawcurr';
dcdFreq = 5e-3; % in ns

%nameList = {'at_neg' 'at_pos'};
%nameList = {'at_neg' 'at_pos' 'gc_neg' 'gc_pos'};
nameList = {'at_neg' 'at_pos' 'gc_neg' 'gc_pos' 'gcm_neg' 'gcm_pos' 'no_pos' 'null_pos'};
ionList = {'K' 'Cl' ''};
%ionList = {''};
currList{1} = {'' '1' '2' '3' '4' '5' '6' '7' '8' '9' '10' 'bp0' 'bp1' 'bp2' 'bp3' 'bp4'};
currList{2} = {'' '1' '2' '3' '4' '5' '6' '7' '8' '9' 'bp0' 'bp1' 'bp2' 'bp3' 'bp4'};
currList{3} = {'' '1' '2' '3' '4' '5' '6' 'R6' '7' '8' '9' '10' '11'...
    '12' '13' '14' '15' '16' '17' 'bp0' 'bp1'};
currList{4} = {'' '1' '2' '3' '4' '5' '6' 'R6' '7' '8' '9' '10' '11'...
    '12' '13' '14' '15' '16' '17' 'bp0' 'bp1' 'bp2'};
currList{5} = {'' '1' '2' '3' '4' '5'};
currList{6} = {'' 'R' '1' 'R1' '2' '3' '4' '5'};
currList{7} = {''};
currList{8} = {''};

out = fopen(outName, 'w');

fprintf(out, 'system\t\tcurrent (pA)\t\tse (pA)\t\ttime (ns)');

for k=1:length(ionList)
    fprintf(out, '\nion %s\n', ionList{k});
    disp(sprintf('\nion %s', ionList{k}));
    for j=1:length(nameList)
        data = [];
        
        currSet = currList{j};
        %currSet = currSet(2:end);
        
        for c=1:length(currSet)
            if ~isempty(currSet{c}) && currSet{c}(1) == 'b'
                fileName = sprintf('%s/%scurr%s_%s_%s.dat', dir, inName, ionList{k}, currSet{c}, nameList{j});
            else 
                fileName = sprintf('%s/%scurr%s_basepair%s_%s.dat', dir, inName, ionList{k}, currSet{c}, nameList{j});
            end
            data = [data; dlmread(fileName, ' ')];
            %if isempty(currSet{c}) || strcmp(currSet{c},'bp0')
            %   d = dlmread(fileName, ' ');
            %   data = [data; d(200:end,:)]; 
            %else
            %   data = [data; dlmread(fileName, ' ')];
            %end
        end

        n = length(data(:,1));
        meanCurr = 1e3*mean(data(:,2));
        stdCurr = 1e3*std(data(:,2));
        errCurr = stdCurr/sqrt(n);
        
        rawFile = sprintf('%s%s_%s_%s.dat', rawPrefix, ionList{k}, nameList{j}, dir);
        %dlmwrite(rawFile, data(:,2), ' ');
        
        fprintf(out, '%s  \t\t%+.1f\t\t%.1f\t\t%.1f\n', nameList{j}, meanCurr, errCurr, n*dcdFreq);
        disp(sprintf('%s  \t\t%+.1f\t\t%.1f\t\t%.1f', nameList{j}, meanCurr, errCurr, n*dcdFreq));
    end
end

fclose(out);
