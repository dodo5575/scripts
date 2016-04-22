#awk '{$1 = printf("%s.log.pos", 1); print}' window_anneal.txt
awk '{$1 = $1.".log.pos"; print}' window_anneal.txt > window_anneal.win
awk '{$1 = $1.".log.pos"; print}' window_middling.txt > window_middling.win
awk '{$1 = $1.".log.pos"; print}' window_raw.txt > window_raw.win
