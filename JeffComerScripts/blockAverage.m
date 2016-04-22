% blockAverage.m
function [tavg xavg xerr] = blockAverage(t, x,block)
% block averages with block count block.

n = length(x);
ind0 = 1:block:n; 
ind1 = [ind0(2:end)-1 n];
ind0 = ind0(1:end-1);
ind1 = ind1(1:end-1);

%ind0 = 1:block:n;
%ind1 = [ind0(2:end)-1 n];
%if ind0(end) == ind1(end)
%  ind0 = ind0(1:end-1);
%  ind1 = ind1(1:end-1);
%end

bn = length(ind0);

tavg = zeros(bn,1);
xavg = zeros(bn,1);
xerr = zeros(bn,1);

for j=1:bn
  bt = t(ind0(j):ind1(j));
  bx = x(ind0(j):ind1(j));
  
  samps = length(bx);
  tavg(j) = mean(bt);
  xavg(j) = mean(bx);
  xerr(j) = std(bx)/sqrt(samps);
end
