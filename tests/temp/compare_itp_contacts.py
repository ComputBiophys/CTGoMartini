#!/usr/bin/env python3
"""
通用 ITP Contacts 比较工具

用于比较两个 GROMACS ITP 文件中的 Go-Martini contacts。
支持不同的原子命名格式和 ID 偏移。

用法:
    python compare_itp_contacts.py <file1> <file2> [options]

示例:
    # 基本比较（自动检测格式）
    python compare_itp_contacts.py file1.itp file2.itp
    
    # 指定 ID 偏移
    python compare_itp_contacts.py file1.itp file2.itp --offset 4
    
    # 指定偏移方向（file1_id + offset = file2_id）
    python compare_itp_contacts.py file1.itp file2.itp --offset 4 --direction forward
    
    # 反向偏移（file1_id - offset = file2_id）
    python compare_itp_contacts.py file1.itp file2.itp --offset 4 --direction reverse
    
    # 自定义输出文件
    python compare_itp_contacts.py file1.itp file2.itp -o comparison_result.txt
"""

import re
import sys
import argparse
from pathlib import Path
from typing import Dict, Tuple, Optional, Callable


def parse_molecule_format(filepath: str, atom_prefix: str = "molecule") -> Dict[Tuple[int, int], Tuple[float, float, float]]:
    """
    解析标准 molecule_X 格式的 ITP 文件
    
    Args:
        filepath: ITP 文件路径
        atom_prefix: 原子名称前缀（默认: molecule）
    
    Returns:
        dict: {(id1, id2): (sigma, epsilon, distance), ...}
    """
    contacts = {}
    pattern = re.compile(
        rf'{atom_prefix}_(\d+)\s+{atom_prefix}_(\d+)\s+1\s+([\d.]+)\s+([\d.]+)\s+;.*?(\d+\.?\d*)'
    )
    
    with open(filepath, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('[') or line.startswith(';'):
                continue
            
            match = pattern.match(line)
            if match:
                id1 = int(match.group(1))
                id2 = int(match.group(2))
                sigma = float(match.group(3))
                epsilon = float(match.group(4))
                distance = float(match.group(5))
                
                key = tuple(sorted((id1, id2)))
                contacts[key] = (sigma, epsilon, distance)
    
    return contacts


def parse_gbp_open_format(filepath: str) -> Dict[Tuple[int, int], Tuple[float, float, float]]:
    """
    解析 gbp_open_X 格式的 ITP 文件（带注释中的ID）
    
    Args:
        filepath: ITP 文件路径
    
    Returns:
        dict: {(id1, id2): (sigma, epsilon, distance), ...}
    """
    contacts = {}
    
    # 尝试匹配带注释的格式
    pattern_with_comment = re.compile(
        r'gbp_open_(\d+)\s+gbp_open_(\d+)\s+1\s+([\d.]+)\s+([\d.]+)\s+;\s+(\d+)\s+(\d+)\s+([\d.]+)'
    )
    
    # 简单格式（无注释）
    pattern_simple = re.compile(
        r'gbp_open_(\d+)\s+gbp_open_(\d+)\s+1\s+([\d.]+)\s+([\d.]+)'
    )
    
    with open(filepath, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('[') or line.startswith(';'):
                continue
            
            # 首先尝试带注释的格式
            match = pattern_with_comment.match(line)
            if match:
                id1 = int(match.group(1))
                id2 = int(match.group(2))
                sigma = float(match.group(3))
                epsilon = float(match.group(4))
                distance = float(match.group(7))
                
                # 验证注释中的ID
                comment_id1 = int(match.group(5))
                comment_id2 = int(match.group(6))
                if id1 != comment_id1 or id2 != comment_id2:
                    print(f"  [行 {line_num}] ID不匹配: {id1}-{id2} vs 注释 {comment_id1}-{comment_id2}")
            else:
                # 尝试简单格式
                match = pattern_simple.match(line)
                if match:
                    id1 = int(match.group(1))
                    id2 = int(match.group(2))
                    sigma = float(match.group(3))
                    epsilon = float(match.group(4))
                    # 估算 distance (sigma * 1.122)
                    distance = sigma * 1.122
                else:
                    continue
            
            key = tuple(sorted((id1, id2)))
            contacts[key] = (sigma, epsilon, distance)
    
    return contacts


def auto_detect_format(filepath: str) -> Tuple[Callable, Optional[int]]:
    """
    自动检测 ITP 文件格式
    
    Args:
        filepath: ITP 文件路径
    
    Returns:
        (parser_function, suggested_offset)
    """
    with open(filepath, 'r') as f:
        content = f.read()
    
    # 检测 gbp_open 格式
    if 'gbp_open_' in content:
        print(f"  检测到格式: gbp_open (带注释ID)")
        return parse_gbp_open_format, None
    
    # 检测 molecule 格式
    molecule_match = re.search(r'molecule_(\d+)', content)
    if molecule_match:
        min_id = int(molecule_match.group(1))
        # 检查是否有其他格式
        if min_id == 1:
            print(f"  检测到格式: molecule (ID从1开始)")
            return lambda fp: parse_molecule_format(fp, "molecule"), 4
        else:
            print(f"  检测到格式: molecule (ID从{min_id}开始)")
            return lambda fp: parse_molecule_format(fp, "molecule"), None
    
    # 默认使用 molecule 格式
    print(f"  未检测到特定格式，使用默认: molecule")
    return lambda fp: parse_molecule_format(fp, "molecule"), None


def shift_contact_ids(contacts: Dict[Tuple[int, int], Tuple], offset: int, direction: str = "forward") -> Dict[Tuple[int, int], Tuple]:
    """
    偏移 contacts 的 ID
    
    Args:
        contacts: 原始 contacts 字典
        offset: 偏移量
        direction: 'forward' (id + offset) 或 'reverse' (id - offset)
    
    Returns:
        偏移后的 contacts 字典
    """
    shifted = {}
    for (id1, id2), value in contacts.items():
        if direction == "forward":
            new_key = tuple(sorted((id1 + offset, id2 + offset)))
        else:
            new_key = tuple(sorted((id1 - offset, id2 - offset)))
        shifted[new_key] = value
    return shifted


def compare_contacts(
    contacts1: Dict[Tuple[int, int], Tuple],
    contacts2: Dict[Tuple[int, int], Tuple],
    file1_name: str,
    file2_name: str,
    tolerance: float = 0.001,
    offset_applied: bool = False
) -> Dict:
    """
    比较两组 contacts
    
    Returns:
        比较结果字典
    """
    keys1 = set(contacts1.keys())
    keys2 = set(contacts2.keys())
    
    common = keys1 & keys2
    only_in_1 = keys1 - keys2
    only_in_2 = keys2 - keys1
    
    # 检查参数差异
    param_diffs = []
    for key in common:
        sigma1, eps1, dist1 = contacts1[key]
        sigma2, eps2, dist2 = contacts2[key]
        
        sigma_diff = abs(sigma1 - sigma2)
        dist_diff = abs(dist1 - dist2)
        
        if sigma_diff > tolerance or dist_diff > tolerance:
            param_diffs.append({
                'key': key,
                'file1': (sigma1, eps1, dist1),
                'file2': (sigma2, eps2, dist2),
                'sigma_diff': sigma_diff,
                'dist_diff': dist_diff
            })
    
    return {
        'total1': len(contacts1),
        'total2': len(contacts2),
        'common': len(common),
        'only_in_1': only_in_1,
        'only_in_2': only_in_2,
        'param_diffs': param_diffs,
        'contacts1': contacts1,
        'contacts2': contacts2
    }


def print_comparison_results(result: Dict, file1_name: str, file2_name: str, offset_info: str = ""):
    """打印比较结果"""
    
    print("\n" + "=" * 80)
    print("CONTACTS 比较结果")
    print("=" * 80)
    
    if offset_info:
        print(f"\nID 映射: {offset_info}")
    
    print(f"\n统计:")
    print(f"  {file1_name}: {result['total1']} 个 contacts")
    print(f"  {file2_name}: {result['total2']} 个 contacts")
    print(f"  共有的 contacts: {result['common']}")
    print(f"  仅在 {file1_name} 中: {len(result['only_in_1'])}")
    print(f"  仅在 {file2_name} 中: {len(result['only_in_2'])}")
    
    # 参数差异
    if result['param_diffs']:
        print(f"\n参数差异 ({len(result['param_diffs'])} 个):")
        for diff in result['param_diffs'][:10]:  # 只显示前10个
            key = diff['key']
            s1, e1, d1 = diff['file1']
            s2, e2, d2 = diff['file2']
            print(f"  {key}:")
            print(f"    {file1_name}: sigma={s1:.6f}, dist={d1:.6f}")
            print(f"    {file2_name}: sigma={s2:.6f}, dist={d2:.6f}")
            print(f"    差异: sigma={diff['sigma_diff']:.6f}, dist={diff['dist_diff']:.6f}")
        if len(result['param_diffs']) > 10:
            print(f"    ... 还有 {len(result['param_diffs']) - 10} 个差异")
    else:
        print("\n✓ 所有共有 contacts 的参数一致")
    
    # 仅在 file1 中的 contacts
    if result['only_in_1']:
        print(f"\n仅在 {file1_name} 中的 contacts ({len(result['only_in_1'])} 个):")
        for key in sorted(result['only_in_1'])[:15]:
            sigma, eps, dist = result['contacts1'][key]
            print(f"  {key}: sigma={sigma:.6f}, eps={eps}, dist={dist:.6f}")
        if len(result['only_in_1']) > 15:
            print(f"  ... 还有 {len(result['only_in_1']) - 15} 个")
    
    # 仅在 file2 中的 contacts
    if result['only_in_2']:
        print(f"\n仅在 {file2_name} 中的 contacts ({len(result['only_in_2'])} 个):")
        for key in sorted(result['only_in_2'])[:15]:
            sigma, eps, dist = result['contacts2'][key]
            print(f"  {key}: sigma={sigma:.6f}, eps={eps}, dist={dist:.6f}")
        if len(result['only_in_2']) > 15:
            print(f"  ... 还有 {len(result['only_in_2']) - 15} 个")
    
    print("\n" + "=" * 80)


def save_detailed_results(result: Dict, file1_name: str, file2_name: str, output_file: str, offset_info: str = ""):
    """保存详细结果到文件"""
    
    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("ITP Contacts 比较详细结果\n")
        f.write("=" * 80 + "\n\n")
        
        f.write(f"文件1: {file1_name}\n")
        f.write(f"文件2: {file2_name}\n")
        if offset_info:
            f.write(f"ID 映射: {offset_info}\n")
        f.write("\n")
        
        f.write("统计:\n")
        f.write(f"  文件1 contacts: {result['total1']}\n")
        f.write(f"  文件2 contacts: {result['total2']}\n")
        f.write(f"  共有 contacts: {result['common']}\n")
        f.write(f"  仅在文件1中: {len(result['only_in_1'])}\n")
        f.write(f"  仅在文件2中: {len(result['only_in_2'])}\n\n")
        
        if result['param_diffs']:
            f.write("-" * 80 + "\n")
            f.write(f"参数差异 ({len(result['param_diffs'])} 个):\n")
            for diff in result['param_diffs']:
                key = diff['key']
                s1, e1, d1 = diff['file1']
                s2, e2, d2 = diff['file2']
                f.write(f"  {key}:\n")
                f.write(f"    文件1: sigma={s1:.6f}, dist={d1:.6f}\n")
                f.write(f"    文件2: sigma={s2:.6f}, dist={d2:.6f}\n")
                f.write(f"    差异: sigma={diff['sigma_diff']:.6f}, dist={diff['dist_diff']:.6f}\n")
            f.write("\n")
        
        if result['only_in_1']:
            f.write("-" * 80 + "\n")
            f.write(f"仅在文件1中的 contacts:\n")
            for key in sorted(result['only_in_1']):
                sigma, eps, dist = result['contacts1'][key]
                f.write(f"  {key}: sigma={sigma:.6f}, eps={eps}, dist={dist:.6f}\n")
            f.write("\n")
        
        if result['only_in_2']:
            f.write("-" * 80 + "\n")
            f.write(f"仅在文件2中的 contacts:\n")
            for key in sorted(result['only_in_2']):
                sigma, eps, dist = result['contacts2'][key]
                f.write(f"  {key}: sigma={sigma:.6f}, eps={eps}, dist={dist:.6f}\n")
            f.write("\n")


def main():
    parser = argparse.ArgumentParser(
        description='比较两个 ITP 文件中的 Go-Martini contacts',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  %(prog)s file1.itp file2.itp
  %(prog)s file1.itp file2.itp --offset 4
  %(prog)s file1.itp file2.itp --offset 4 --direction reverse
  %(prog)s file1.itp file2.itp -o result.txt
        """
    )
    
    parser.add_argument('file1', help='第一个 ITP 文件')
    parser.add_argument('file2', help='第二个 ITP 文件')
    parser.add_argument('--offset', type=int, default=None,
                        help='ID 偏移量（默认：自动检测）')
    parser.add_argument('--direction', choices=['forward', 'reverse'], default='forward',
                        help='偏移方向：forward (file1_id + offset = file2_id) 或 reverse (file1_id - offset = file2_id)')
    parser.add_argument('--tolerance', type=float, default=0.001,
                        help='参数比较的容差（默认：0.001）')
    parser.add_argument('-o', '--output', default='comparison_result.txt',
                        help='输出文件路径（默认：comparison_result.txt）')
    parser.add_argument('--format1', choices=['molecule', 'gbp_open', 'auto'], default='auto',
                        help='文件1的格式（默认：自动检测）')
    parser.add_argument('--format2', choices=['molecule', 'gbp_open', 'auto'], default='auto',
                        help='文件2的格式（默认：自动检测）')
    
    args = parser.parse_args()
    
    # 检查文件存在
    for f in [args.file1, args.file2]:
        if not Path(f).exists():
            print(f"错误: 文件不存在: {f}")
            sys.exit(1)
    
    file1_path = Path(args.file1)
    file2_path = Path(args.file2)
    
    print("=" * 80)
    print("ITP Contacts 比较工具")
    print("=" * 80)
    print(f"\n文件1: {file1_path}")
    print(f"文件2: {file2_path}")
    
    # 解析文件1
    print(f"\n解析 {file1_path.name}...")
    if args.format1 == 'auto':
        parser1, suggested_offset1 = auto_detect_format(args.file1)
        contacts1 = parser1(args.file1)
    elif args.format1 == 'molecule':
        contacts1 = parse_molecule_format(args.file1)
        suggested_offset1 = 4
    else:
        contacts1 = parse_gbp_open_format(args.file1)
        suggested_offset1 = None
    print(f"  找到 {len(contacts1)} 个 contacts")
    
    # 解析文件2
    print(f"\n解析 {file2_path.name}...")
    if args.format2 == 'auto':
        parser2, suggested_offset2 = auto_detect_format(args.file2)
        contacts2 = parser2(args.file2)
    elif args.format2 == 'molecule':
        contacts2 = parse_molecule_format(args.file2)
    else:
        contacts2 = parse_gbp_open_format(args.file2)
    print(f"  找到 {len(contacts2)} 个 contacts")
    
    if not contacts1 or not contacts2:
        print("\n错误: 未能解析到 contacts")
        sys.exit(1)
    
    # 确定偏移量
    offset = args.offset
    if offset is None:
        # 尝试自动确定偏移量
        min_id1 = min(min(k) for k in contacts1.keys())
        min_id2 = min(min(k) for k in contacts2.keys())
        offset = min_id2 - min_id1
        print(f"\n自动计算 ID 偏移: {offset} (基于最小ID: {min_id1} -> {min_id2})")
    else:
        print(f"\n使用指定 ID 偏移: {offset}")
    
    # 应用偏移
    offset_info = ""
    if offset != 0:
        if args.direction == 'forward':
            contacts1_shifted = shift_contact_ids(contacts1, offset, 'forward')
            offset_info = f"{file1_path.name}_id + {offset} = {file2_path.name}_id"
        else:
            contacts1_shifted = shift_contact_ids(contacts1, offset, 'reverse')
            offset_info = f"{file1_path.name}_id - {offset} = {file2_path.name}_id"
        print(f"偏移方向: {args.direction}")
    else:
        contacts1_shifted = contacts1
        offset_info = "无偏移"
    
    # 比较
    print(f"\n正在比较...")
    result = compare_contacts(
        contacts1_shifted, contacts2,
        file1_path.name, file2_path.name,
        tolerance=args.tolerance
    )
    
    # 打印结果
    print_comparison_results(result, file1_path.name, file2_path.name, offset_info)
    
    # 保存结果
    output_path = Path(args.output)
    save_detailed_results(result, file1_path.name, file2_path.name, str(output_path), offset_info)
    print(f"详细结果已保存到: {output_path.absolute()}")
    
    # 返回退出码（如果有差异则返回1）
    if result['only_in_1'] or result['only_in_2'] or result['param_diffs']:
        print("\n注意: 发现差异")
        return 1
    else:
        print("\n✓ 两个文件完全一致")
        return 0


if __name__ == "__main__":
    sys.exit(main())
