
# source double-precision gmx_d
source /usr/local/gromacs-2020.7-d/bin/GMXRC
set -e
# use double precision for everything
gmx_d grompp -f martini_md.mdp -c /home/ys/CommonUse/Martini/test/CTGoMartini/ctgomartini/tests/api/../data/MBGoMartini/Beta2AR/Strfiles/Beta2AR_No9.gro -p system.top -o Beta2AR_No9.tpr -r /home/ys/CommonUse/Martini/test/CTGoMartini/ctgomartini/tests/api/../data/MBGoMartini/Beta2AR/Strfiles/Beta2AR_No9.gro -maxwarn 1 -v
gmx_d mdrun -deffnm Beta2AR_No9 -rerun /home/ys/CommonUse/Martini/test/CTGoMartini/ctgomartini/tests/api/../data/MBGoMartini/Beta2AR/Strfiles/Beta2AR_No9.gro -nt 1 -v
echo pot | gmx_d energy -f Beta2AR_No9 -o energy_Beta2AR_No9
echo 0 | gmx_d traj -f Beta2AR_No9 -s Beta2AR_No9 -of forces_Beta2AR_No9
rm mdout.mdp
rm Beta2AR_No9.edr Beta2AR_No9.log Beta2AR_No9.tpr Beta2AR_No9.trr
