function ret = bondDistribRandom(n, r0, spring, kT)

rc = 0.5*r0 + 0.5*sqrt(r0^2 + 8*kT/spring);
pc = 4*pi*rc^2*exp(-0.5*spring*(rc-r0)^2/kT);

rMean = rc;
rStd = sqrt(kT/(0.8*spring));

ret = zeros(1,n);
j = 0;
while j<n
    r = rMean + rStd*randn;
    s = pc*exp(-0.5*((r-rMean)/rStd)^2);
    p = 4*pi*r^2*exp(-0.5*spring*(r-r0)^2/kT);
    
    if r > 0 && p/s > rand
        j = j + 1;
        ret(j) = r;
    end
end
