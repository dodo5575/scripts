clear all;

%list1 = [649.8		13.7		373.6];
%list2 = [643.0		13.9		370.7];

list1 = [620.7		8.3		1048.2];
list2 = [+614.6		9.3		863.5];

t1 = list1(3);
di1 = list1(2);
i1 = list1(1);

t2 = list2(3);
di2 = list2(2);
i2 = list2(1);

curr = (t1*i1 + t2*i2)/(t1+t2)
err = sqrt((t1/(t1+t2)*di1)^2 + (t2/(t1+t2)*di2)^2)
tim = t1 + t2
