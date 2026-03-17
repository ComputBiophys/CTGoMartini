
# 修正 resid 和 segid
echo "修正 resid 和 segid 信息..."
python ChangeResidSegid_input.py -f npt.gro -s "A:1-271;A:1-271" -o npt_revised.pdb || echo "处理 npt.pdb 失败"
python ChangeResidSegid_input.py -f Open/Open_cg.pdb -s "A:1-271;A:1-271" -o Open/Open_cg.pdb || echo "处理 Up_cg.pdb 失败"
python ChangeResidSegid_input.py -f Closed/Closed_cg.pdb -s "A:1-271;A:1-271" -o Closed/Closed_cg.pdb || echo "处理 Down_cg.pdb 失败"

# dRMS calculation
python dRMS_analysis_nc.py -s npt_revised.pdb -nc output.nc --ref_str Open/Open_cg.pdb Closed/Closed_cg.pdb --num_workers 10 -o dRMStraj_nc --chunk 4
