function [Ave,Error] = blockAverage(data, blockSize)
    
Nlength = length(data);

reducedData = [];

if mod(Nlength,blockSize) ~= 0
    
    Nblock = floor(Nlength / blockSize);
    
    for i = 0:(Nblock-1)
        
        reducedData = [reducedData;mean(data(i * blockSize + 1: i * blockSize + blockSize))];
    
    end
    
    reducedData = [reducedData;mean(data(Nblock * blockSize + 1: end))];
    
    
else
    
    Nblock = Nlength / blockSize;
    
    for i = 0:(Nblock-1)
        
        reducedData = [reducedData;mean(data(i * blockSize + 1: i * blockSize + blockSize))];
    
    end
    
    
end

meanText = sprintf('The average is %d', mean(reducedData));
%disp(meanText)

Ave = mean(reducedData);

stdText = sprintf('The standard deviation is %d', std(reducedData));
%disp(stdText)

Error = std(reducedData);

