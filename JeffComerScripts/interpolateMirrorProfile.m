clear all;

% Parameters:
surfList = {'raw' 'middling' 'anneal'};
profileFile = sprintf('diffuse_dev_dmmp25_DoubleCut.0.dat');
surfZ0List = [23.65970588 24.30987952 23.30121457];
surfZ1List = [55.65029411 58.61012048 58.91878542];
templateZ0 = 23.65970588;
waterDiffuse = 530;
molDiffuse = 91;
outPrefix = 'dmmp_diffuse_';

plotColor = [0 0 0; 0.9 0 0; 0.8 0.4 0; 0.5 0.5 0; 0 0.8 0; 0 0.7 0.7; 0 0 0.9; 0.6 0 0.6; 0.6 0.6 0.6];

figure(1)
clf
hold on

for j=1:length(surfList)
    surf = surfList{j};
    
    % Load the template profile.
    templateFile = sprintf('md_%s.z.dat', surf);    
    template = dlmread(templateFile, ' ');
    templateZ = template(:,1);
    templateV = template(:,2);

    % Load the profile.
    data = dlmread(profileFile, ' ');
    z = data(:,1);
    v = data(:,2);
    
    % Map the values.
    templateSurf = pchip(z, v, templateZ);
    n = length(templateZ);
    for k=1:n
        if templateZ(k) > z(end)
            templateSurf(k) = waterDiffuse;
        elseif templateZ(k) < z(1)
            templateSurf(k) = v(1);
        end
    end
    
    % Mirror.
    mirrorZ = 0.5*(surfZ0List(j) + surfZ1List(j));
    dz = templateZ(2) - templateZ(1);
    mirrorInd = round((mirrorZ - templateZ(1))/dz) + 1;
    for k=mirrorInd:n
        if n - k < 1
            templateSurf(k) = v(1); 
        else
            templateSurf(k) = templateSurf(n - k);
        end
    end
    % Scale to the molecule's diffusivity.
    templateSurf = templateSurf*molDiffuse/waterDiffuse;
    
    % Write the diffusion profile.
    outFile = sprintf('%s%s.dat', outPrefix, surf);
    dlmwrite(outFile, [templateZ templateSurf], ' ');
    
    % Plot.
    color0 = plotColor(mod(j-1,length(plotColor))+1,:);
    gh = plot(templateZ, templateSurf);
    set(gh, 'Color', color0);
    set(gh, 'Marker', 'none');
    set(gh, 'MarkerSize', 7);
    set(gh, 'LineStyle', '-');
end

hold off