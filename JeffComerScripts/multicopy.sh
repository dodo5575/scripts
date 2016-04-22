for file in $(find . -iname "*1400*.curr"); do
/bin/cp -a $file /projects/jcomer/brownian/diffusion_map/spool
echo $file
done
