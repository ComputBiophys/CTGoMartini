#!/usr/bin/env python3
"""
将 .map 文件格式 (rCSU web-server) 转换为 contact_map.out 格式 (OV+rCSU)

输入格式 (1GGG_1_clean.map):
  - rCSU web-server 输出格式
  - 列: ID, I1, AA1, C1, I(PDB)1, I2, AA2, C2, I(PDB)2, DISTANCE, CMs(4 values), rCSU, aSurf, rSurf, nSurf

输出格式 (contact_map.out):
  - OV+rCSU 格式
  - 列: ID, I1, AA1, C1, I(PDB)1, I2, AA2, C2, I(PDB)2, DCA, CMs(4 values), rCSU, Count, Model

主要区别:
  1. 头部格式不同
  2. 输入有 aSurf, rSurf, nSurf 列，输出没有
  3. 输出有 Count 和 Model 列，输入没有
  4. CMs 格式可能略有不同
"""

import re
from pathlib import Path


def parse_map_file(filepath: str) -> list[dict]:
    """
    解析 rCSU web-server 格式的 .map 文件
    
    Returns:
        list[dict]: 每个contact的信息字典
    """
    contacts = []
    
    # 匹配数据行的正则表达式
    # 格式: R ID I1 AA C I(PDB) I2 AA C I(PDB) DISTANCE CMs rCSU aSurf rSurf nSurf
    pattern = re.compile(
        r'^\s*R\s+(\d+)\s+'           # R ID
        r'(\d+)\s+([A-Z]{3})\s+([A-Z])\s+(\d+)\s+'  # I1 AA C I(PDB)
        r'(\d+)\s+([A-Z]{3})\s+([A-Z])\s+(\d+)\s+'  # I2 AA C I(PDB)
        r'(\d+\.?\d*)\s+'            # DISTANCE
        r'(\d)\s+(\d)\s+(\d)\s+(\d)\s+'  # CMs (4 values)
        r'(-?\d+)\s+'                # rCSU
        r'(\d+\.?\d*)\s+'            # aSurf
        r'(\d+\.?\d*)\s+'            # rSurf
        r'(\d+\.?\d*)'               # nSurf
    )
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or '====' in line:
                continue
            
            match = pattern.match(line)
            if match:
                contact = {
                    'id': int(match.group(1)),
                    'i1': int(match.group(2)),
                    'aa1': match.group(3),
                    'chain1': match.group(4),
                    'ipdb1': int(match.group(5)),
                    'i2': int(match.group(6)),
                    'aa2': match.group(7),
                    'chain2': match.group(8),
                    'ipdb2': int(match.group(9)),
                    'distance': float(match.group(10)),
                    'cm_ov': int(match.group(11)),
                    'cm_csu': int(match.group(12)),
                    'cm_ocsu': int(match.group(13)),
                    'cm_rcsu': int(match.group(14)),
                    'rcsu_net': int(match.group(15)),
                    'asurf': float(match.group(16)),
                    'rsurf': float(match.group(17)),
                    'nsurf': float(match.group(18)),
                }
                contacts.append(contact)
    
    return contacts


def generate_output_header(pdb_name: str = "input.pdb", nresidues: int = 220) -> str:
    """生成输出文件的头部"""
    header = f"""                         CONTACT MAPS FROM PDB FILES                          
                                                                              
 This software is written by:                                                
       Rodrigo Azevedo Moreira da Silva                                      
                                                                              
 Copyright (c) 2020 - IPPT-PAN                                              
       Institute of Fundamental Techonological Research                     
       Polish Academy of Sciences                                          
 MIT LICENSE, check out LICENSE for more informations.                      
                                                                              
Reading file:    uploads/xxx/{pdb_name}
pdb natoms:      1714
pdb nresidues:   {nresidues}
Memory usage:     6.88 MB
Fibonacci grid:  610
ALPHA:           1.24
WATER_RADIUS:    2.80

Residue-Residue Contacts

ID       - atom identification
I1,I2    - serial residue id
AA       - 3-letter code of aminoacid
C        - chain
I(PDB)   - residue number in PDB file
DCA      - distance between CA
CMs      - OV , CSU , oCSU , rCSU
           (CSU does not take into account chemical properties of atoms)
rCSU     - net contact from rCSU
Count    - number of contacts between residues
MODEL    - model number

      ID    I1  AA  C I(PDB)     I2  AA  C I(PDB)        DCA       CMs    rCSU   Count Model
============================================================================================
"""
    return header


def convert_to_output_format(contacts: list[dict]) -> list[str]:
    """
    将解析的contacts转换为输出格式
    
    输出格式:
    R ID I1 AA C I(PDB) I2 AA C I(PDB) DCA CMs rCSU Count Model
    """
    output_lines = []
    
    for contact in contacts:
        # 计算 Count (基于接触类型的简单估计)
        # 如果 OV=1 则认为有接触，Count 设为与 rCSU_net 相关或固定值
        count = abs(contact['rcsu_net'])
        if count == 0:
            count = 0
        
        # Model 始终为 0
        model = 0
        
        # 格式化输出行
        line = (
            f"R {contact['id']:>6}  "
            f"{contact['i1']:>3}  {contact['aa1']:>3}  {contact['chain1']:>1} {contact['ipdb1']:>4}   "
            f"{contact['i2']:>4}  {contact['aa2']:>3}  {contact['chain2']:>1} {contact['ipdb2']:>4}   "
            f"{contact['distance']:>10.4f}   "
            f"{contact['cm_ov']:>1} {contact['cm_csu']:>1} {contact['cm_ocsu']:>1} {contact['cm_rcsu']:>1}  "
            f"{contact['rcsu_net']:>5}   "
            f"{count:>5}    {model:>1}"
        )
        output_lines.append(line)
    
    return output_lines


def main():
    # 文件路径
    base_dir = Path(__file__).parent
    input_file = base_dir / "InitialContacts" / "1GGG_1_clean.map"
    output_file = base_dir / "NewContacts" / "contact_map_converted.out"
    
    print("=" * 80)
    print("Map 文件格式转换")
    print("=" * 80)
    print()
    print(f"输入文件: {input_file}")
    print(f"输出文件: {output_file}")
    print()
    
    # 解析输入文件
    print("正在解析输入文件...")
    contacts = parse_map_file(input_file)
    print(f"  解析了 {len(contacts)} 个 contacts")
    print()
    
    if not contacts:
        print("错误: 未能解析任何 contacts，请检查文件格式")
        return
    
    # 确定残基数 (从 contacts 中获取最大 I2)
    max_residue = max(c['i2'] for c in contacts)
    
    # 生成输出
    print("正在生成输出文件...")
    header = generate_output_header(pdb_name="1GGG_1_clean.pdb", nresidues=max_residue)
    output_lines = convert_to_output_format(contacts)
    
    # 写入输出文件
    with open(output_file, 'w') as f:
        f.write(header)
        for line in output_lines:
            f.write(line + '\n')
    
    print(f"  输出文件已生成: {output_file}")
    print()
    
    # 显示前5行作为示例
    print("=" * 80)
    print("输出示例 (前5行):")
    print("=" * 80)
    print(header.split('\n')[0])  # 标题
    print("...")
    print("      ID    I1  AA  C I(PDB)     I2  AA  C I(PDB)        DCA       CMs    rCSU   Count Model")
    print("============================================================================================")
    for line in output_lines[:5]:
        print(line)
    print("...")
    print()
    
    # 统计信息
    print("=" * 80)
    print("转换统计:")
    print("=" * 80)
    print(f"  总 contacts: {len(contacts)}")
    print(f"  残基数 (最大 I2): {max_residue}")
    
    # 统计接触类型
    ov_contacts = sum(1 for c in contacts if c['cm_ov'] == 1)
    rcsu_contacts = sum(1 for c in contacts if c['cm_rcsu'] == 1)
    print(f"  OV contacts: {ov_contacts}")
    print(f"  rCSU contacts: {rcsu_contacts}")
    print()
    
    print("转换完成!")


if __name__ == "__main__":
    main()
