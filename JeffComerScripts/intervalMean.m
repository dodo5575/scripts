function [meanVal seVal nVal] = intervalMean(ti, tf, t, val)
% [meanVal seVal nVal] = intervalMean(ti, tf, t, val)
% Returns the mean, standard error, and number of points of data
% within the t-interval ti to tf. The data has ordinates t and 
% values val.
% Author: Jeff Comer <jcomer2@illinois.edu>

nVal = 0;
sumVal = 0;
sumValSq = 0;
for step=1:length(t)
    % Find all points in this time interval.
    if ti <= t(step) && tf > t(step)
        sumVal = sumVal + val(step);
        sumValSq = sumValSq + val(step)*val(step);
        nVal = nVal + 1;
    end
end

if nVal == 0
    meanVal = 0;
    seVal = 0;
    return
end

meanVal = sumVal/nVal;
meanValSq = sumValSq/nVal;
seVal = sqrt(meanValSq - meanVal*meanVal)/nVal;



