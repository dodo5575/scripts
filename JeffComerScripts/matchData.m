% matchData.m
clear all;

%file0 = 'number/num_bam_specific_1.5V.dat';
%file1 = 'electrode/top_bam_specific_1.5V.dat';

%file0 = 'number/num_bam_spec_nogrid_2V.dat';
%file1 = 'electrode/top_bam_spec_nogrid_2V.dat';

%file0 = 'number/num_eco_start_4V.dat';
%file1 = 'electrode/top_eco_start_4V.dat';

%file0 = 'number/num_eco_start_2.5V.dat';
%file1 = 'electrode/top_eco_start_2.5V.dat';

%file0 = 'number/num_eco_grisha_2Vb.dat';
%file1 = 'electrode/top_eco_grisha_2Vb.dat';

%file0 = 'number/num_eco_grisha_orig_2V_100ps.dat';
%file1 = 'electrode/top_eco_grisha_orig_2V_100ps.dat';

%file0 = 'number/num_bam_specific_1.5V.dat';
%file1 = 'electrode/top_bam_specific_1.5V.dat';

%file0 = 'number/num_bam_specific_2V.dat';
%file1 = 'electrode/top_bam_specific_2V.dat';

%file0 = 'number/num_bam_nonspec_1.5V.dat';
%file1 = 'electrode/top_bam_nonspec_1.5V.dat';

%file0 = 'number/num_bam_nonspec_1Va.dat';
%file1 = 'electrode/top_bam_nonspec_1Va.dat';

%file0 = 'number/num_bam_nonspec_0.5V.dat';
%file1 = 'electrode/top_bam_nonspec_0.5V.dat';

data0 = dlmread(file0, ' ');
data1 = dlmread(file1, ' ');

start = 1;
data1 = data1(start:end,:);

dt = 1e-6*1.0*5000*10;
t = data1(:,1) - 0.5*dt;
g = data1(:,2);
n = length(t);
f = zeros(n,1);
n0 = length(data0(:,1));

for j=1:n
    k = binsearch(t(j), data0(:,1));
    if k >= n0
       f(j) = data0(end,2);
       continue
    end
    if k <= 0
       f(j) = data0(1,2);
       continue
    end
    
    t0 = data0(k,1);
    t1 = data0(k+1,1);
    if abs(t(j)-t0) <= abs(t(j)-t1)
        f(j) = data0(k,2);
    else
        f(j) = data0(k+1,2);
    end
end

N = length(f);
sumF = sum(f);
sumFF = sum(f.*f);
sumG = sum(g);
sumFG = sum(f.*g);

file0
alpha = (sumFG - sumF*sumG/N)/(sumFF - sumF*sumF/N)
beta = (sumG - alpha*sumF)/N

fNew = alpha*f + beta;

%figure(1)
%plot(t, f)
%figure(2)
%plot(t, g)
figure(3)
plot(t, g, 'ko', t, fNew, 'b-')



