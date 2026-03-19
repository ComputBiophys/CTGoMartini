
# source double-precision gmx_d
source /usr/local/gromacs-2020.7-d/bin/GMXRC
set -e
# use double precision for everything
gmx_d grompp -f martini_md.mdp -c /home/ys/CommonUse/Martini/test/CTGoMartini/ctgomartini/tests/api/../data/MBGoMartini/GlnBP/Strfiles/GlnBP_No9.gro -p system.top -o GlnBP_No9.tpr -r /home/ys/CommonUse/Martini/test/CTGoMartini/ctgomartini/tests/api/../data/MBGoMartini/GlnBP/Strfiles/GlnBP_No9.gro -maxwarn 1 -v
gmx_d mdrun -deffnm GlnBP_No9 -rerun /home/ys/CommonUse/Martini/test/CTGoMartini/ctgomartini/tests/api/../data/MBGoMartini/GlnBP/Strfiles/GlnBP_No9.gro -nt 1 -v
echo pot | gmx_d energy -f GlnBP_No9 -o energy_GlnBP_No9
echo 0 | gmx_d traj -f GlnBP_No9 -s GlnBP_No9 -of forces_GlnBP_No9
rm mdout.mdp
rm GlnBP_No9.edr GlnBP_No9.log GlnBP_No9.tpr GlnBP_No9.trr
