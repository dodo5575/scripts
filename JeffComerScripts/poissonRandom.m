function k = poissonRandom(n0)
l = exp(-n0);

k = 0;
p = rand;
while p >= l
    p = p * rand;
    k = k + 1;
end
