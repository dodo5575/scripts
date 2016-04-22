function [a b da db] = linearRegression(x, y)
% Get the linear regression.
% b is the slope
% a is the offset
% Author: Jeff Comer <jcomer2@illinois.edu>% 

% Compute the sums.
sum1 = length(x);
sumX = sum(x);
sumY = sum(y);
sumXX = sum(x.*x);
sumXY = sum(x.*y);

% Do the linear regression.
delta = sum1*sumXX - sumX*sumX;
a = (sumXX*sumY - sumX*sumXY)/delta;
b = (sum1*sumXY - sumX*sumY)/delta;

% Compute the uncertainties.
dy = sqrt(1/(sum1-2)*sum((y - a - b*x).^2));
da = dy*sqrt(sumXX/delta);
db = dy*sqrt(sum1/delta);
