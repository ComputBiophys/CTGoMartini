# Note:
# MDTraj time: There is a problem in which the time series of trjactory are not the same as the input one.

python run_ctgomartini.py -i npt.inp
python run_REMD.py -i md.inp

python Distance_nc4.py -s npt_revised.pdb --skip 1 --num_workers 20

python Distance_nc.py -s npt_revised.pdb --skip 1 --num_workers 20


python ExtractStateItp.py -f system.top -m GlnBP -s 1 -o GlnBP_stateA.itp
python ExtractStateItp.py -f system.top -m GlnBP -s 2 -o GlnBP_stateB.itp

cp system.top system_stateA.top
sed -i 's/GlnBP.itp/GlnBP_stateA.itp/g' system_stateA.top

cp system.top system_stateB.top
sed -i 's/GlnBP.itp/GlnBP_stateB.itp/g' system_stateB.top

python run_REMD.py


python ChangeResidSegid_input.py -f npt.pdb -s "A:5-224;A:5-224" -o npt_revised.pdb || echo "处理 npt.pdb 失败"
#python ChangeResidSegid_input.py -f Open/Open_cg.pdb -s "A:5-224;A:5-224" -o Open/Open_cg.pdb || echo "处理 Open_cg.pdb 失败"
#python ChangeResidSegid_input.py -f Closed/Closed_cg.pdb -s "A:5-224;A:5-224" -o Closed/Closed_cg.pdb || echo "处理 Closed_cg.pdb 失败"

python Distance_nc.py -s npt_revised.pdb -f "replica_*.xtc" -o distance_nc.dat -p 10
