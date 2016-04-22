clear all;

fprintf('\n****%s\n', datestr(now));

curr0 = 5.95456750479;
err0 = 0.0659425147025;
curr1 = 2.34266221913;
err1 = 0.0658434511332;

frac = curr1/curr0
err = abs(1.0/curr0)*err1 + abs(-curr1/curr0^2)*err0

curr0 = 5.95456750479;
err0 = 0.0659425147025;
curr1 = 2.12920445587;
err1 = 0.0912516660947;

frac = curr1/curr0
err = abs(1.0/curr0)*err1 + abs(-curr1/curr0^2)*err0
