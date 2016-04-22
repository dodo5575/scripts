for x in anneal anneal_exclude middling middling_exclude raw1 raw1_exclude phant phant_exclude; do
	./condorSubmission.pl $x.brown 4
	condor_submit $x.brown.sub
done
