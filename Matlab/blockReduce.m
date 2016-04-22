function reducedData = blockReduce(data, blockSize)
    
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


