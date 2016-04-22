clear all;

outFile = 'other_current.dat';
outTable = 'other_table.dat';
prefix = 'thru90_curr_';
suffix = '.dat';
inList = {'2V' [5 -1] '3V' [5 -1] '4V' [5 -1]...
    'open_enzyme_pore_2V' [4 -1]...
    'force_dna_1V' [0 -1] 'free_dna_0V' [0 -1] 'free_dna_1V' [0 -1] 'free_dna_2V' [8 -1]...
    'eco_grisha_2Vb' [15 -1] 'eco_grisha_orig_2V_100ps' [0 40] 'eco_start_2.5V' [8 12] 'eco_start_4V' [8 -1]...
    'bam_nonspec_0.5V' [8 -1] 'bam_nonspec_1Va' [8 -1] 'bam_nonspec_1.5V' [8 -1] 'bam_2Va' [5 -1] ...
    'bam_specific_1V' [10 -1] 'bam_specific_1.5V' [20 30] 'bam_specific_2V' [5 12]...
    'bam_spec_nogrid_2V' [8 20] 'bam_spec_grid_2V' [5 -1]...
    'open_enzyme_hot_2V' [4 -1] 'dna_mem10_1V' [8 -1] 'dna_mem5_1V' [8 -1] 'str_dna_1V' [8 18]...
    'open_trap2.0_1V' [3 -1] 'trap2.0_1V' [3 -1]...
    'period_trap_0Va' [10 -1] 'period_trap_0.5V' [10 -1] 'period_trap_0.5V' [15 -1] 'period_trap_1Va' [8 22]};

%inList = {'2V' [5 -1]};

out = fopen(outFile, 'w');
outTable = fopen(outTable, 'w');

inNum = length(inList)/2;
for j=0:inNum-1
    name = inList{2*j + 1};
    cut = inList{2*j + 2};
    
    fileName = sprintf('%s%s%s', prefix, name, suffix);
    
    data = dlmread(fileName, ' ');
    t = data(:,1);
    curr = data(:,2);
    
    % Find the places to cut the data.
    cut0 = -1;
    cut1 = length(t);
    for s=1:length(t)
        if cut0 < 0 && t(s) >= cut(1)
            cut0 = s;
        end
        
        if cut(2) > 0 && t(s) < cut(2)
            cut1 = s;
        end
    end
    if cut0 < 1
        cut0 = 1;
    end
   
    mat = regexp(name, '[0123456789\.]+V', 'match');
    mat1 = mat{1};
    voltage = str2double(mat1(1:end-1));
    
    t = t(cut0:cut1);
    curr = curr(cut0:cut1);
    n = length(t);
    
    meanTime = mean(t);
    meanCurr = mean(curr);
    meanCond = meanCurr/voltage;
    stdCond = std(curr/voltage);
    errCond = stdCond/sqrt(n);
    stdCurr = std(curr);
    errCurr = stdCurr/sqrt(n);
    totalTime = t(end) - t(1);

    if regexp(name, 'eco')
        enzyme = '\eco';
        seq = 'cognate';
    elseif regexp(name, 'bam_spec')
        enzyme = '\bam';
        seq = 'cognate';
    elseif regexp(name, 'bam')
        enzyme = '\bam';
        seq = 'nonspec.';
    elseif regexp(name, 'force') 
        enzyme = 'none';
        seq = 'fixed';
    elseif regexp(name, 'free')
        enzyme = 'none';
        seq = 'free';
    elseif regexp(name, 'dna')
        enzyme = 'none';
        seq = 'fixed';
    else
        enzyme = 'none';
        seq = 'none';
    end
    
    % Write the table file.
    fprintf(outTable, '%s & %s & %.1f & %0.2f $\\pm$ %0.2f & %0.2f $\\pm$ %0.2f \\\\\n', ...
        enzyme, seq, voltage, meanCurr, errCurr, meanCond, errCond);
    
    % Write the text file.
    fprintf(out, 'name: %s\n', name);
    fprintf(out, 'cutTime0: %g\n', cut(1));
    fprintf(out, 'cutTime1: %g\n', cut(2));
    fprintf(out, 'meanTime: %g\n', meanTime);
    fprintf(out, 'totalTime: %g\n', totalTime);
    fprintf(out, 'meanCurrent: %g\n', meanCurr);
    fprintf(out, 'errCurrent: %g\n', errCurr);
    fprintf(out, 'voltage: %g\n', voltage);
    fprintf(out, 'meanConductance: %g\n', meanCond);
    fprintf(out, 'errConductance: %g\n', errCond);
    fprintf(out, '\n', errCond);
end
fclose(out);
fclose(outTable);
