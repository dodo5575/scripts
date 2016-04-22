clear all;

% Parameters:
run = 'three';
coulombConst = 566.44320919705/92;
charge = [1 -1];
radius = [1.76375 2.27];
eps = [0.0870 0.15];
nameList = {'pot-pot' 'chl-chl' 'pot-chl'};
selList = [1 1; 2 2; 1 2];
dir = '.';
%z = linspace(2,11,200);

for j=1:length(nameList) 
    s1 = selList(j,1);
    s2 = selList(j,2);
    name = nameList{j};
    
    inName = sprintf('%s/pmf_%s_%s.dx.z.dat', dir, run, name);
    data = dlmread(inName, ' ');
    
    r = data(:,1);
    v = data(:,2);
    u = coulombConst*charge(s1)*charge(s2)./r + computeHardcore(r,eps(s1), eps(s2), radius(s1), radius(s2));
    uc = coulombConst*charge(s1)*charge(s2)./r;
    
    cut0 = floor(length(r)/4);
    cut1 = 3*cut0;
    meanV = mean(v(cut0:cut1));
    meanU = mean(u(cut0:cut1));
    v1 = v - meanV + meanU;
    
    outNameU = sprintf('%s/%s_energy_%s.dat', dir, run, name);
    outNameUc = sprintf('%s/%s_coulomb_%s.dat', dir, run, name);
    outNameV = sprintf('%s/%s_free_%s.dat', dir, run, name);
        
    dlmwrite(outNameU, [r u], ' ');
    fprintf('Wrote %s.\n', outNameU);
    dlmwrite(outNameUc, [r uc], ' ');
    fprintf('Wrote %s.\n', outNameUc);
    dlmwrite(outNameV, [r v1], ' ');
    fprintf('Wrote %s.\n', outNameV);
end
