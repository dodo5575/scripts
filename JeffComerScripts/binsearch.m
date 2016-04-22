function a = binsearch(value, sortedList)
% Perform a binary search to find the indices of the values that
% immediately surround VALUE in sorted list SORTEDLIST.
% For SORTEDLIST(1) <= VALUE < SORTEDLIST(end),
% the result is A is such that SORTEDLIST(A) <= VALUE < SORTEDLIST(A+1).
% If VALUE is less than all members of the list, 0 is returned.
% If VALUE is greater than all members of the list, the last index is
% returned.
% The complexity is log(n).
% Author: Jeff Comer <jcomer2@illinois.edu>

low = 1;
high = length(sortedList);

% Check if the value falls outside the list.
% Note that the return values do not contain valid indices.
if value < sortedList(low)
    a = 0;
    return
end
if value >= sortedList(high)
    a = high;
    return
end

while (high - low > 1)
    mid = floor((high+low)/2);
    
    if (value < sortedList(mid))
        high = mid;
    else
        low = mid;
    end
    
    %[low high]
end

a = low;
return



