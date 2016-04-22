#!/bin/bash

bibCmd='{10}\\setlength{\\itemsep}{-1mm}'
prefix=part
section=0
for name in significance objectives; do
    # Make the bibfiles for this section.
    f=${prefix}_${name}
    pdflatex $f
    bibtex $f

    pre="S${section}_"

    # Change names of bibitems using sed.
    # Also make the items closer together.
    sed -e s/\bibitem{/\bibitem{${pre}/ -e s/{10}/${bibCmd}/ -e s/{1}/${bibCmd}/ -e s/Http/http/ -e s/")\."/")"/ < $f.bbl > ${f}_secbib.bbl

    # Refactor the citations using the prefix.
	    tclsh prefixLatexCite.tcl $name.tex $pre ${name}_secbib.tex
	    
    # Increment the section number.
	    ((section++))
done

# We have to do it differently for the computation section, since it's
# broken into parts.
for name in computation; do
    # Make the bibfiles for this section.
    f=${prefix}_${name}
    pdflatex $f
    bibtex $f

    pre="S${section}_"

    # Change names of bibitems using sed.
    sed -e s/\bibitem{/\bibitem{${pre}/ -e s/{10}/${bibCmd}/ -e s/{1}/${bibCmd}/ -e s/Http/http/ -e s/")\."/")"/ < $f.bbl > ${f}_secbib.bbl

    # Refactor the citations using the prefix.
	    tclsh prefixLatexCite.tcl $name.tex $pre ${name}_secbib.tex
	    for sf in approach job-characterization parallelPerformance development workflow io data; do
		tclsh prefixLatexCite.tcl $sf.tex $pre ${sf}_secbib.tex
	    done
	    
    # Increment the section number.
	    ((section++))
done
