index=$1
gmx grompp -f em.mdp -c window_${index}.pdb -r window_${index}.pdb -p ../../system.top -o em${index}.tpr
gmx mdrun -v -deffnm em${index} -nt 8

gmx grompp -f npt.mdp -c em${index}.gro -r window_${index}.pdb -p ../../system.top -o npt${index}.tpr -n umbrella.ndx
gmx mdrun -v -deffnm npt${index} -nt 8

gmx grompp -f md.mdp -c npt${index}.gro -r window_${index}.pdb -p ../../system.top -o md${index}.tpr -n umbrella.ndx
gmx mdrun -v -deffnm md${index} -nt 8

