#!/usr/bin/env python3
"""
比较两个ITP文件中的contacts (版本2)

文件1: NewContacts2/go_nbparams.itp - 使用 molecule_X 格式, ID从1开始
文件2: InitialContacts/gbp_open_go-table_VirtGoSites.itp - 使用 gbp_open_X 格式, ID从5开始

ID映射关系: 文件2的ID = 文件1的ID + 4
"""

import re
from pathlib import Path


def parse_go_nbparams(filepath: str) -> dict:
    """
    解析 go_nbparams.itp 文件
    格式: molecule_X molecule_Y 1 sigma epsilon ;go bond distance
    
    Returns:
        dict: {(id1, id2): (sigma, epsilon, distance), ...}
    """
    contacts = {}
    pattern = re.compile(
        r'molecule_(\d+)\s+molecule_(\d+)\s+1\s+([\d.]+)\s+([\d.]+)\s+;go bond\s+([\d.]+)'
    )
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('['):
                continue
            
            match = pattern.match(line)
            if match:
                id1 = int(match.group(1))
                id2 = int(match.group(2))
                sigma = float(match.group(3))
                epsilon = float(match.group(4))
                distance = float(match.group(5))
                
                # 标准化contact键（小id在前，大id在后）
                key = tuple(sorted((id1, id2)))
                contacts[key] = (sigma, epsilon, distance)
    
    return contacts


def parse_virtgosites(filepath: str) -> dict:
    """
    解析 gbp_open_go-table_VirtGoSites.itp 文件
    格式: gbp_open_X gbp_open_Y 1 sigma epsilon ; id1 id2 distance
    
    Returns:
        dict: {(id1, id2): (sigma, epsilon, distance), ...}
    """
    contacts = {}
    # 首先尝试匹配带注释的格式
    pattern_with_comment = re.compile(
        r'gbp_open_(\d+)\s+gbp_open_(\d+)\s+1\s+([\d.]+)\s+([\d.]+)\s+;\s+(\d+)\s+(\d+)\s+([\d.]+)'
    )
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(';'):
                continue
            
            match = pattern_with_comment.match(line)
            if match:
                id1 = int(match.group(1))
                id2 = int(match.group(2))
                sigma = float(match.group(3))
                epsilon = float(match.group(4))
                # 分号后的id和distance是注释信息
                comment_id1 = int(match.group(5))
                comment_id2 = int(match.group(6))
                distance = float(match.group(7))
                
                # 验证注释中的id是否与前面的id一致
                if id1 != comment_id1 or id2 != comment_id2:
                    print(f"警告: 行中的ID不匹配: {line[:80]}...")
                
                # 标准化contact键（小id在前，大id在后）
                key = tuple(sorted((id1, id2)))
                contacts[key] = (sigma, epsilon, distance)
    
    return contacts


def shift_contact_ids(contacts: dict, offset: int) -> dict:
    """
    将所有contact的ID进行偏移
    
    Args:
        contacts: {(id1, id2): value, ...}
        offset: 偏移量（新id = 原id + offset）
    
    Returns:
        dict: 偏移后的contacts
    """
    shifted = {}
    for (id1, id2), value in contacts.items():
        new_key = tuple(sorted((id1 + offset, id2 + offset)))
        shifted[new_key] = value
    return shifted


def compare_contacts(raw_contacts: dict, initial_contacts: dict, tolerance: float = 0.001):
    """
    比较两组contacts
    
    Args:
        raw_contacts: 从go_nbparams.itp解析的contacts
        initial_contacts: 从gbp_open_go-table_VirtGoSites.itp解析的contacts
        tolerance: 距离比较的容差
    """
    # 将raw_contacts的ID偏移+4，以便与initial_contacts比较
    raw_shifted = shift_contact_ids(raw_contacts, offset=4)
    
    raw_keys = set(raw_shifted.keys())
    initial_keys = set(initial_contacts.keys())
    
    # 分类contacts
    common = raw_keys & initial_keys  # 共有的
    only_in_raw = raw_keys - initial_keys  # 只在raw中
    only_in_initial = initial_keys - raw_keys  # 只在initial中
    
    print("=" * 80)
    print("CONTACTS 比较结果")
    print("=" * 80)
    print()
    print(f"NewContacts2/go_nbparams.itp (ID偏移+4后): {len(raw_shifted)} 个contacts")
    print(f"InitialContacts/gbp_open_go-table_VirtGoSites.itp: {len(initial_contacts)} 个contacts")
    print()
    
    # 统计
    print("-" * 80)
    print("统计:")
    print(f"  共有的 contacts: {len(common)}")
    print(f"  仅在 NewContacts2 中: {len(only_in_raw)}")
    print(f"  仅在 InitialContacts 中: {len(only_in_initial)}")
    print()
    
    # 检查共有contacts的参数差异
    if common:
        print("-" * 80)
        print("共有 contacts 的参数差异 (显示sigma/distance差异超过容差的):")
        diff_count = 0
        for key in sorted(common):
            raw_sigma, raw_epsilon, raw_dist = raw_shifted[key]
            init_sigma, init_epsilon, init_dist = initial_contacts[key]
            
            sigma_diff = abs(raw_sigma - init_sigma)
            dist_diff = abs(raw_dist - init_dist)
            
            if sigma_diff > tolerance or dist_diff > tolerance:
                diff_count += 1
                print(f"  Contact {key}:")
                print(f"    NewContacts2:    sigma={raw_sigma:.6f}, dist={raw_dist:.6f}, eps={raw_epsilon}")
                print(f"    InitialContacts: sigma={init_sigma:.6f}, dist={init_dist:.6f}, eps={init_epsilon}")
                print(f"    差异: sigma_diff={sigma_diff:.6f}, dist_diff={dist_diff:.6f}")
        
        if diff_count == 0:
            print("  所有共有contacts的参数在容差范围内一致")
        else:
            print(f"  共有 {diff_count} 个contacts存在参数差异")
        print()
    
    # 仅在NewContacts2中的contacts
    if only_in_raw:
        print("-" * 80)
        print(f"仅在 NewContacts2/go_nbparams.itp 中的 contacts ({len(only_in_raw)} 个):")
        for key in sorted(only_in_raw)[:20]:  # 只显示前20个
            sigma, epsilon, dist = raw_shifted[key]
            print(f"  {key}: sigma={sigma:.6f}, dist={dist:.6f}")
        if len(only_in_raw) > 20:
            print(f"  ... 还有 {len(only_in_raw) - 20} 个")
        print()
    
    # 仅在InitialContacts中的contacts
    if only_in_initial:
        print("-" * 80)
        print(f"仅在 InitialContacts/gbp_open_go-table_VirtGoSites.itp 中的 contacts ({len(only_in_initial)} 个):")
        for key in sorted(only_in_initial)[:20]:  # 只显示前20个
            sigma, epsilon, dist = initial_contacts[key]
            print(f"  {key}: sigma={sigma:.6f}, dist={dist:.6f}")
        if len(only_in_initial) > 20:
            print(f"  ... 还有 {len(only_in_initial) - 20} 个")
        print()
    
    # 保存详细结果到文件
    output_file = Path(__file__).parent / "contact_comparison_result_v2.txt"
    with open(output_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("CONTACTS 比较详细结果 (NewContacts2 vs InitialContacts)\n")
        f.write("=" * 80 + "\n\n")
        
        f.write("文件:\n")
        f.write("  文件1: NewContacts2/go_nbparams.itp\n")
        f.write("  文件2: InitialContacts/gbp_open_go-table_VirtGoSites.itp\n")
        f.write("  ID映射: 文件1的ID + 4 = 文件2的ID\n\n")
        
        f.write("统计:\n")
        f.write(f"  文件1 contacts数: {len(raw_shifted)}\n")
        f.write(f"  文件2 contacts数: {len(initial_contacts)}\n")
        f.write(f"  共有 contacts: {len(common)}\n")
        f.write(f"  仅在文件1中: {len(only_in_raw)}\n")
        f.write(f"  仅在文件2中: {len(only_in_initial)}\n\n")
        
        if only_in_raw:
            f.write("-" * 80 + "\n")
            f.write(f"仅在 NewContacts2/go_nbparams.itp 中的 contacts:\n")
            for key in sorted(only_in_raw):
                sigma, epsilon, dist = raw_shifted[key]
                f.write(f"  {key}: sigma={sigma:.6f}, dist={dist:.6f}\n")
            f.write("\n")
        
        if only_in_initial:
            f.write("-" * 80 + "\n")
            f.write(f"仅在 InitialContacts/gbp_open_go-table_VirtGoSites.itp 中的 contacts:\n")
            for key in sorted(only_in_initial):
                sigma, epsilon, dist = initial_contacts[key]
                f.write(f"  {key}: sigma={sigma:.6f}, dist={dist:.6f}\n")
            f.write("\n")
    
    print("=" * 80)
    print(f"详细结果已保存到: {output_file}")
    print("=" * 80)


def main():
    # 文件路径
    base_dir = Path(__file__).parent
    file1 = base_dir / "NewContacts2" / "go_nbparams.itp"
    file2 = base_dir / "InitialContacts" / "gbp_open_go-table_VirtGoSites.itp"
    
    print("=" * 80)
    print("比较两个ITP文件的contacts (NewContacts2 vs InitialContacts)")
    print("=" * 80)
    print()
    print(f"文件1: {file1}")
    print(f"文件2: {file2}")
    print()
    
    # 解析两个文件
    print("正在解析文件...")
    raw_contacts = parse_go_nbparams(file1)
    print(f"  从 NewContacts2/go_nbparams.itp 解析了 {len(raw_contacts)} 个contacts")
    
    initial_contacts = parse_virtgosites(file2)
    print(f"  从 gbp_open_go-table_VirtGoSites.itp 解析了 {len(initial_contacts)} 个contacts")
    print()
    
    # 比较
    compare_contacts(raw_contacts, initial_contacts)


if __name__ == "__main__":
    main()
