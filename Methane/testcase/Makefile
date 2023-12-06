%.tpr: ../%.gro
	gmx grompp -f init.mdp -c $< -p CH4.top -o $*.tpr -maxwarn 9999999 > $*.grompp.log
%-last.gro: %.tpr
	gmx mdrun -s $*.tpr -c $*-last.gro -deffnm $* > $*.mdrun.log
