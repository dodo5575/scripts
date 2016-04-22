clear all;

% Parameters:
run = 'hydration';
dir = '.';
rMax = 30;
points = 400;
coulombConst = 332.06369; % A/e^2 kcal/mol
dielectric = 92; % TIP3P water
kT = 0.5862292; % kcal/mol
charge = [1 -1];
nameList = {'pot-chl' 'pot-pot' 'chl-chl'};
selList = [1 2; 1 1; 2 2];

% Lennard-Jones paramters from D. Beglovd and B. Roux (1994).
radius = [1.76375 2.27];
eps = [0.0870 0.15];
% hydration potential parameters from W. Im and B. Roux (2002).
c0List = [-3.7 -0.60 -0.50];
c1List = [2.90 4.40 4.90];
c2List = [0.90 0.90 0.90];
c3List = [0.75 0.80 0.080];
c4List = [0 0.25 0.25];

r = linspace(0, rMax, points);
plotColor = [0 0 0; 0.9 0 0; 0.8 0.4 0; 0.5 0.5 0; 0 0.8 0; 0 0.7 ...
	     0.7; 0 0 0.9; 0.6 0 0.6; 0.6 0.6 0.6];

figure(1)
clf
hold on

for j=1:length(nameList) 
    s1 = selList(j,1);
    s2 = selList(j,2);
    name = nameList{j};
    
    potElec = coulombConst/dielectric*charge(s1)*charge(s2)./r;
    potLj = computeLennardJones(r, eps(s1), eps(s2), radius(s1), radius(s2));
    potHydr = computeHydrationPotential(r, c0List(j), c1List(j), c2List(j), c3List(j), c4List(j));
    pot = (potElec + potLj + potHydr)/kT;
    
    color0 = plotColor(mod(j-1,length(plotColor))+1,:);
    gh = plot(r,pot);
    set(gh, 'Color', color0);
    set(gh, 'Marker', 'o');
    set(gh, 'MarkerSize', 7);
    set(gh, 'LineStyle', '-');
    
    outName = sprintf('%s/hydration_%s.dat', dir, name);
    dlmwrite(outName, [r pot], ' ');
end

hold off
axis([0 rMax/2 -4 4])
