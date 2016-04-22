rm log-local-err.log log-local-out.log
qsub -q tlong -S /bin/sh -pe namd 48 -e log-local-err.log -o log-local-out.log -cwd -N loop runLoop.sh
