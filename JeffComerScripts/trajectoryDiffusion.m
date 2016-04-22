clear all;

inDir = 'traj_chloride';
prefix = 'follow_basepair_no_1M_neg*';
cutIndex = 40;
dimensions = 2;

figure(1)
clf
hold on

glob = sprintf('%s/%s.dat', inDir, prefix);
fileList = dir(glob);
nFiles = length(fileList);
fprintf('Found %d files.\n', nFiles);

varSum = zeros(cutIndex,1);
varCount = 0;
for j=1:nFiles
    fileName = sprintf('%s/%s', inDir, fileList(j).name);
    trajName = sprintf('%s/xy_%s', inDir, fileList(j).name);
    trajName1 = sprintf('%s/trsq_%s', inDir, fileList(j).name);
    
    data = dlmread(fileName, ' ');
    
    if length(data(:,1)) < cutIndex
        continue
    end

    dt = data(2,1) - data(1,1);
    n = length(data(:,1));
    t = dt*((1:n)-1)';
    %z = data(:,2);
    %z = z - z(1);
    r = data(:,2:4);
    
    %r = [r(:,1)-r(1,1) r(:,2)-r(1,2) r(:,3)-r(1,3)];
    %Ignore z since we've biased the distribution in z by rejecting
    % data that leaves a certain region.
    r = [r(:,1)-r(1,1) r(:,2)-r(1,2)];
    rsq = dot(r,r,2);
    
    varSum = varSum + rsq(1:cutIndex);
    varCount = varCount + 1;
    
    x = r(:,1);
    y = r(:,2);
    dlmwrite(trajName, [x y], ' ');
    dlmwrite(trajName1, [t rsq], ' ');

    %plot(t,rsq, 'k.-')
end
hold off

varTime = ((1:cutIndex)'-1)*dt;
varMean = varSum/varCount;
[a b da db] = linearRegression(varTime, varMean);
varMeanFit = a + b*varTime;
diffusion = b/(2*dimensions);
diffusionErr = db/(2*dimensions);

fprintf('diffusion const: %.2f +/- %.1f\n', diffusion, diffusionErr);

figure(2)
plot(varTime, varMean, 'bo-', varTime, varMeanFit, 'k--');
xlabel('t (ns)');
ylabel('<x^2 + y^2> A^2');
dispText = sprintf('slope = %.2f, D = %.2f +/- %.2f', b, diffusion, diffusionErr);
text(0.02, 160, dispText);
