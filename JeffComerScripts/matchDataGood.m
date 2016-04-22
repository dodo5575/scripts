% matchData.m
clear all;

name = {'eco_start_4V' 'eco_start_2.5V' 'eco_grisha_2Vb' 'eco_grisha_orig_2V_100ps'...
    'bam_specific_2V' 'bam_spec_nogrid_2V' 'bam_specific_1.5V' 'bam_nonspec_1.5V' 'bam_nonspec_1Va'...
    'bam_nonspec_0.5V'};

dir0 = 'dna_stretch';
prefix0 = 'str_';
dir1 = 'electrode_good';
prefix1 = 'top_';
%len0 = 0.34;

for j=9:9%1:length(name)
    file0 = sprintf('%s/%s%s.dat', dir0, prefix0, name{j});
    file1 = sprintf('%s/%s%s.dat', dir1, prefix1, name{j});

    data0 = dlmread(file0, ' ');
    data1 = dlmread(file1, ' ');

    start = 1;
    data0 = data0(start:end,:);
    data1 = data1(start:end,:);

    t = data1(:,1);
    f = 1./(data0(:,2)/10);
    g = data1(:,2);

    N = length(f);
    sumF = sum(f);
    sumFF = sum(f.*f);
    sumG = sum(g);
    sumFG = sum(f.*g);

    file0
    alpha = (sumFG - sumF*sumG/N)/(sumFF - sumF*sumF/N)
    beta = (sumG - alpha*sumF)/N
    %alpha = -7.4
    beta = 0

    %fNew = alpha*f + beta;
    gNew = (g-beta)/alpha;

    figure(4)
    plot(t, f, 'k-', t, gNew, 'r-')
    %plot(t, f, 'k-')
    
    outF = sprintf('invlen_%s.dat', name{j});
    outG = sprintf('transpot_%s.dat', name{j});
    outTrans = sprintf('trans_%s.dat', name{j});
    dlmwrite(outF, [t f], ' ');
    dlmwrite(outG, [t gNew], ' ');
    dlmwrite(outTrans, [alpha; beta], ' ');
end



