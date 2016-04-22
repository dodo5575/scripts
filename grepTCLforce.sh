#1 (kcal/mol angstrom) = 69.43521595 pN

#This script extracts the TCLforce output from log file.
#Usage: ./grepTCLforce.sh name logFile startTime dimension
#By Chen-Yu Li	cli56@illinois.edu
#2013/7/16

case "$4" in 

    x)  
	grep "calculated force of sampleDNA" $2 | awk "{print \$2*0.000002 + $3, \$8*69.43521595}" > ${1}_X_force.dat
	;;

    y)  
	grep "calculated force of sampleDNA" $2 | awk "{print \$2*0.000002 + $3, \$9*69.43521595}" > ${1}_Y_force.dat
	;;

    z)  
	grep "calculated force of sampleDNA" $2 | awk "{print \$2*0.000002 + $3, \$10*69.43521595}" > ${1}_Z_force.dat
	;;

esac





