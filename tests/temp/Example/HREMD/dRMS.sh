
# 修正 resid 和 segid
echo "修正 resid 和 segid 信息..."
python ChangeResidSegid_input.py -f npt.pdb -s "A:40-315;B:40-315;A:40-315;B:40-315" -o npt_revised.pdb || echo "处理 npt.pdb 失败"
python ChangeResidSegid_input.py -f Up/Up_cg.pdb -s "A:40-315;B:40-315;A:40-315;B:40-315" -o Up/Up_cg.pdb || echo "处理 Up_cg.pdb 失败"
python ChangeResidSegid_input.py -f Down/Down_cg.pdb -s "A:40-315;B:40-315;A:40-315;B:40-315" -o Down/Down_cg.pdb || echo "处理 Down_cg.pdb 失败"

# 查找所有副本轨迹文件（支持多位数编号）
echo "查找副本轨迹文件..."
replica_files=(replica_[0-9]*.xtc)  # 匹配 replica_1.xtc, replica_11.xtc, ...
if [ ${#replica_files[@]} -eq 0 ]; then
    echo "在 $dir 中未找到副本文件"
    cd ..
    continue
fi

# 处理每个副本文件
for replica_file in "${replica_files[@]}"; do
    # 提取纯数字编号（去掉 "replica_" 和 ".xtc"）
    j=${replica_file#replica_}  # 去掉前缀 → "1.xtc", "11.xtc"
    j=${j%.xtc}                # 去掉后缀 → "1", "11"

    echo "处理 $replica_file (副本 $j)"

    # 计算 dRMS
    python MBAnalysis_v4_argparse.py -s npt_revised.pdb -f "$replica_file" \
        -r Up/Up_cg.pdb Down/Down_cg.pdb -sel 'name BB' \
        -prefix "dRMStraj_replica_${j}" -n 20 || echo "处理 $replica_file 失败"
done

