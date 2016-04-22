% meanOnInterval.m
function [v dv] = meanOnInterval(x, x0, x1, y)

n = length(x);

count = 0;
sumY = 0.0;
sumYY = 0.0;
for j=1:n
    if x(j) >= x0 && x(j) < x1
        count = count + 1;
        sumY = sumY + y(j);
	sumYY = sumYY + y(j)*y(j);
    end
end

if count == 0
    v = 0;
    dv = 0;
else
    v = sumY/count;
    dv = sqrt((sumYY-sumY*sumY/count)/(count-1)/count);
end
