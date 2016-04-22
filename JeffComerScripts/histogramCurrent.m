% histogramCurrent.m
clear all;

%name = 'open_p1.6_4V0-6';
%name = 'open_p2.0_4V0-2';
%name = 'open_p1.2_4V0-3';
%name = 'eloop_p2.0_4V0-0';
name = 'ecoil_first_4V0-2';
%name = 'open_4V1_4V-21';
type = '';

inName = sprintf('curr%s_%s.dat', type, name);
outName = sprintf('histo%s_%s.dat', type, name);
outNameProb = sprintf('prob%s_%s.dat', type, name);
outNameMean = sprintf('mean%s_%s.dat', type, name);

cut = 400;

data = dlmread(inName, ' ');
t = data(cut:end,1);
curr = data(cut:end,2);

meanCurr = mean(curr)
stdCurr = std(curr)
seCurr = std(curr)/sqrt(length(curr))

[n x] = hist(curr, 14);
area = sum(n)*(x(2)-x(1));
prob = n/area;


plot(x, prob)
dlmwrite(outName, [x' n'], ' ');
dlmwrite(outNameProb, [x' prob'], ' ');
dlmwrite(outNameMean, [meanCurr seCurr stdCurr], ' ');



