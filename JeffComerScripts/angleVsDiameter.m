clear all;

name = {'coil_p1.2_6.5Va0-1' 'coil_p1.2_4Va0-5' 'coilx_p1.6_6.5Va0-2'...
    'coilx_p1.6_6.5V0-3' 'coilx_p1.6_4V0-14' 'coil_p2.0_2V0-7' 'coil_p2.0_4V0-2'};
suff = '.20.block';

plotColor = [0 0 0; 0.9 0 0; 0.8 0.4 0; 0.5 0.5 0; 0 0.8 0; 0 0.7 0.7;...
    0 0 0.9; 0.5 0.3 0.8; 0.6 0 0.6; 0.5 0.5 0.5];
figure(1)
clf
hold on

disp(sprintf('\nangleVsDiameter'))
for j=1:length(name)
    posFile = sprintf('last1_pos_%s.dat%s', name{j}, suff);
    diamFile = sprintf('diam_%s.dat%s', name{j}, suff);
    angleFile = sprintf('last1_base_angle_%s.dat%s', name{j}, suff);
    outFilePos = sprintf('angleVpos_%s.dat', name{j});
    outFileDiam = sprintf('angleVdiam_%s.dat', name{j});
    
    % Read the data.
    data = dlmread(posFile, ' ');
    z = data(:,2);
    data = dlmread(diamFile, ' ');
    d = data(:,2);
    data = dlmread(angleFile, ' ');
    phi = data(:,2);
    
    % Write the angle v. position.
    dlmwrite(outFilePos, [z phi], ' ');
    
    % Remove the negative positions.
    n = 0;
    dPlus = zeros(size(z));
    phiPlus = zeros(size(z));
    for k=1:length(z)
        if z(k) >= 0.0
            n = n + 1;
            dPlus(n) = d(k);
            phiPlus(n) = phi(k);
        end
    end
    dPlus = dPlus(1:n);
    phiPlus = phiPlus(1:n);
    disp(sprintf('Processed: %s %d', name{j}, length(phiPlus)));
    
    % Write the angle v. diameter.
    dlmwrite(outFileDiam, [dPlus phiPlus], ' ');
    
    % Plot the data.
    color0 = plotColor(mod(j-1,length(plotColor))+1,:);
    %gh = plot(z, phi);
    gh = plot(dPlus, phiPlus);
    set(gh, 'Color', color0);
    set(gh, 'Marker', 'o');
    set(gh, 'MarkerSize', 10);
    set(gh, 'LineStyle', '-');
end
hold off
h = legend(name);
set(h, 'Interpreter', 'none');




