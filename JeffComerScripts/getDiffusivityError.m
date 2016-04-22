clear all;

inDir = '.';

rankList = {'piece0' 'piece1' 'total'};
totalRank = 3;
regionList = {'pot_all' 'chl_all'};
outFile = 'error_and_diffusion.txt';
sliceList = [0];
slicePosList = [0];
sliceWidth = 0.15;
cutTime = 0.02;
dimensions = 2;
%convFactor = 10; % from A^2/ns to um^2/s
convFactor = 1.0; % results A^2/ns


plotColor = [0 0 0; 0.9 0 0; 0.8 0.4 0; 0.5 0.5 0; 0 0.8 0; 0 0.7 0.7; 0 0 0.9; 0.6 0 0.6; 0.6 0.6 0.6];
figure(1)
clf
hold on
count = 0;

out = fopen(outFile, 'w');
fprintf('\n****%s\n', datestr(now));
%fprintf(out, '\n****%s\n', datestr(now));


for reg=1:length(regionList)
    region = regionList{reg};
    for slice=1:length(sliceList)
        sliceNum = sliceList(slice);
        slicePos = slicePosList(slice);
        
        diffusion = zeros(length(rankList),1);
        slope = zeros(size(diffusion));
        inter = zeros(size(diffusion));
        msd = cell(length(rankList), 1);
        
        for ra=1:length(rankList)
            rank = rankList{ra};
            
            % Read the data.
            fileName = sprintf('%s_dev_%s.0.%d.msd', rank, region, sliceNum);
            data = dlmread(fileName, ' ');
            fprintf('%s\n', fileName);
            
            if isempty(data)
                continue
            end
            dt = data(2,1) - data(1,1);
            cutIndex = floor(cutTime/dt);
            if cutIndex < 1
                cutIndex = 1;
            end
            
            if length(data(:,1)) < cutIndex
                continue
            end
          
            % Get the data.
            t = data(cutIndex:end,1);
            n = length(t);
            msd{ra} = data(cutIndex:end,2);       
            
            % Compute the linear regression.
            [a b da db] = linearRegression(t, msd{ra});
            slope(ra) = b;
            inter(ra) = a;
            
            % Compute the diffusivity.
            diffusion(ra) = b/(2*dimensions);            
        end
        
        % Compute the theoretical msd using the total.
        theo = inter(totalRank) + slope(totalRank)*t;
        diffusionTotal = diffusion(totalRank);
        totalMsd = msd{totalRank};
        
        % Add the legend entries.
        name = sprintf('%s_%d', region, sliceNum);
        count = count + 1;
        legendList{2*(count-1)+1} = name;
        legendList{2*(count-1)+2} = name;
        
        % Compute the best diffusion and error and convert.
        diffTotal = convFactor*diffusionTotal;
        diffErr = convFactor*range(diffusion);
        diffErr0 = diffTotal - convFactor*min(diffusion);
        diffErr1 = convFactor*max(diffusion) - diffTotal;
        
        fprintf('%s: %.3f -%.3f +%.3f\n', name, diffTotal, diffErr0, diffErr1);
        fprintf(out, '%s: %.3f -%.3f +%.3f\n', name, diffTotal, diffErr0, diffErr1);
        
        %fprintf(out, '%.3g %.3g %.3g %.3g\n', slicePos, diffTotal, diffErr0, diffErr1);
        
        % Plot it.
        color0 = plotColor(mod(count-1,length(plotColor))+1,:);
        gh = plot(t, totalMsd);
        set(gh, 'Color', color0);
        set(gh, 'Marker', '.');
        set(gh, 'MarkerSize', 7);
        set(gh, 'LineStyle', 'None');
        
        gh = plot(t, theo);
        set(gh, 'Color', color0);
        set(gh, 'LineStyle', '-');
        
        % Say the result.
        %dispText = sprintf('slope = %.3f, D = %.3f +/- %.3f', b, diffusion, diffErr);
        %text(0.1, 4+0.5*j, dispText);
        % Output in um^2/s
    end
end

fclose(out);

lh = legend(legendList);
set(lh, 'Interpreter', 'none');
xlabel('t (ns)');
ylabel('<x^2 + y^2> A^2');
hold off

