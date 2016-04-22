function result = gaussianFilter(t, x, sigma)
% Read a .dx file.
% Y = GAUSSIANFILTER(T,X,SIGMA) applies a gaussian filter to the data
% X with a standard deviation SIGMA.
% Author: Jeff Comer <jcomer2@illinois.edu>

sigmaSq = sigma*sigma;
n = length(t);
dt = (t(end) - t(1))/n;
disp(sprintf('Average data interval: %.3g', dt));

window = ceil(8.0*sigma/dt);
windowHalf = ceil(window/2);
disp(sprintf('Using a window of four standard deviations: %.3g', window));

result = zeros(size(t));
index = 1;
for j=1:n
    start = index - windowHalf;
    stop = index + windowHalf;
    
    if start < 1
        start = 1;
    end
    
    if stop > n
        stop = n;
    end
    
    w = exp(-(t(j)-t(start:stop)).*(t(j)-t(start:stop))/(2.0*sigmaSq));
    sumWeights = sum(w);
    result(j) = sum(x(start:stop).*w)/sumWeights;

    index = index + 1;
end



