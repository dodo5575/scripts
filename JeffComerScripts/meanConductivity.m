% meanCurrent.m
clear all;

outName ='cg';
prefix = 'cg_bulk_';
inName = {'0.03' '0.1' '0.2' '0.4' '0.56' '0.8' '1.0'};
%outName ='martini';
%prefix = 'martini_bulk_';
%inName = {'0.03' '0.1' '0.2' '0.56' '1.0'};

suffix = '_1V.dat';
cutTime = 3; % starting time (ns)
dcdFreq = 1000;
timeStep = 40.0; % fs
outFile = sprintf('curr_conc_%s.dat', outName);
outFile1 = sprintf('conduct_conc_%s.dat', outName);
outFile2 = sprintf('lambda_conc_%s.dat', outName);
outFile3 = sprintf('fit_conc_%s.dat', outName);

l = 81.035; % angstroms
voltage = 1; % volts
lambda0 = 7.5549;

n = length(inName);
out = zeros(n, 3);
out1 = zeros(n, 3);
out2 = zeros(n, 3);
out3 = zeros(n, 3);

x = zeros(n, 1);
y = zeros(n, 1);
for j=1:n 
  %dataK = dlmread(sprintf('currK_%s',inName), ' ');
  %dataCl = dlmread(sprintf('currCl_%s',inName), ' ');
  fileName = sprintf('curr_%s%sM%s', prefix, inName{j}, suffix);
  data = dlmread(fileName, ' ');
  conc = str2double(inName{j});
  
  step = timeStep*dcdFreq*1e-6;
  cut = ceil(cutTime/step);

  t = data(cut:end,1);
  curr = data(cut:end,2);
  samples = length(t);
  
  meanCurr = mean(curr);
  stdCurr = std(curr);
  errCurr = std(curr)/sqrt(samples);
  
  cond = curr*1e-9/(voltage*l*1e-10);
  meanCond = mean(cond);
  errCond = std(cond)/sqrt(samples);
  
  lambda = cond./conc;
  meanLambda = mean(lambda);
  errLambda = std(lambda)/sqrt(samples);
  
  disp(sprintf('\nCurrent for concentration of %s M:', inName{j}))
  disp(sprintf('mean: %.10g', meanCurr))
  disp(sprintf('std: %.10g', stdCurr)) 
  disp(sprintf('err: %.10g', errCurr))
  disp(sprintf('mean cond: %.10g', meanCond))
  disp(sprintf('err cond: %.10g', errCond))
  
  out(j,1) = conc;
  out(j,2) = meanCurr;
  out(j,3) = errCurr;
  
  out1(j,1) = conc;
  out1(j,2) = meanCond;
  out1(j,3) = errCond;
    
  out2(j,1) = conc;
  out2(j,2) = meanLambda; % mS*m^2/mol
  out2(j,3) = errLambda; % mS*m^/mol
    
  x(j) = conc;
  y(j) = meanLambda;
end

 % Fit.
  %B = (sum(y) - sum(y.*sqrt(x)) * sum(sqrt(x)))/sum(x)/(n-sum(sqrt(x))^2/sum(x));
  %A = (n*B - sum(y))/sum(sqrt(x));
  %yFit = B - A*sqrt(x);
%plot(x, y, 'k-o', x, yFit, 'b-')
errorbar(out2(:,1), out2(:,2), out2(:,3))

dlmwrite(outFile, out, ' ');
dlmwrite(outFile1, out1, ' ');
dlmwrite(outFile2, out2, ' ');
%dlmwrite(outFile3, out3, ' ');



