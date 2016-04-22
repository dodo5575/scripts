#This script extracts the SMD output from log file.
#Usage: ./grepSMDforce.sh name logFile startTime dimension
#By Chen-Yu Li	cli56@illinois.edu
#2013/7/16

case "$4" in 

    x)  
	grep "SMD  " $2 | awk "{print \$2*0.000002 + $3, \$3}" > ${1}_X_position.dat
	grep "SMD  " $2 | awk "{print \$2*0.000002 + $3, \$6}" > ${1}_X_force.dat
	;;

    y)  
	grep "SMD  " $2 | awk "{print \$2*0.000002 + $3, \$4}" > ${1}_Y_position.dat
	grep "SMD  " $2 | awk "{print \$2*0.000002 + $3, \$7}" > ${1}_Y_force.dat
	;;

    z)  
	grep "SMD  " $2 | awk "{print \$2*0.000002 + $3, \$5}" > ${1}_Z_position.dat
	grep "SMD  " $2 | awk "{print \$2*0.000002 + $3, \$8}" > ${1}_Z_force.dat
	;;

esac



