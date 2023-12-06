%.gro: %.cif
	python3 cif2gro.py < $< > $@

clean:
	rm "#"* *.pdb

%-0.gro: %.tpr Makefile
	prev=`expr $* - 1`; prev=`printf '%05d' $$prev`; echo $$prev; \
	if [ -e $$prev.part*.trr ] ; \
		then echo 0 | gmx trjconv -pbc nojump -f $$prev.part*.trr -s $*.tpr -o $*-0.gro; \
	elif [ -e $*.trr ] ; \
		then echo 0 | gmx trjconv -pbc nojump -f $*.trr -s $*.tpr -o $*-0.gro; \
	fi
