clear all;

dcdFreq = 20*5e-3;
cutTime = 3.5;
sys = {'b-dna' 's-dna' 'open'};
reg = {'pore' 'far'};
ion = {'K' 'Cl'};
concFactor = 55.523;

fprintf('\n****%s\n', datestr(now));
fprintf('sys\t\tconc (mol/kg)\t\t err (mol/kg)\n');
cutInd = ceil(cutTime/dcdFreq);

for i=1:length(ion)
    for r=1:length(reg)
        for s=1:length(sys)
            inFile = sprintf('num%s_pore6_%s_500mV_%s.dat', ion{i}, sys{s}, reg{r});
            data = dlmread(inFile, ' ');
            numIon = data(cutInd:end,2);
            
            inFile = sprintf('numWater_pore6_%s_500mV_%s.dat', sys{s}, reg{r});
            data = dlmread(inFile, ' ');
            numWater = data(cutInd:end,2);
            
            nw = mean(numWater);
            ew = std(numWater)/sqrt(length(numWater));
            ni = mean(numIon);
            ei = std(numIon/sqrt(length(numIon)));
            
            concIon = concFactor*ni/nw;
            errIon = concFactor*(1/nw*ei + ni/nw^2*ew);
                        
            fprintf('%s %s %s\t\t', sys{s}, reg{r}, ion{i});
            fprintf('%.3g\t\t%.3g\n', concIon, errIon);
        end
    end
end
