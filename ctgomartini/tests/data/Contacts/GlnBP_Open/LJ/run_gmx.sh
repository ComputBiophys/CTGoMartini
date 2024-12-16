
# source double-precision gmx_d
source /usr/local/gromacs-2020.7-d/bin/GMXRC
set -e
# use double precision for everything
gmx_d grompp -f martini_md.mdp -c minimized.gro -p system.top -o conf.tpr -r minimized.gro -maxwarn 1 -v
gmx_d mdrun -deffnm conf -rerun minimized.gro -nt 1 -v
echo pot | gmx_d energy -f conf -o energy_conf
echo 0 | gmx_d traj -f conf -s conf -of forces_conf
rm mdout.mdp
rm conf.edr conf.log conf.tpr conf.trr
