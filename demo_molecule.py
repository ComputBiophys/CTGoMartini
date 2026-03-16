#!/usr/bin/env python3
"""
Molecule 类使用示例
展示如何创建、填充和初始化分子拓扑
"""

from ctgomartini.core import Molecule, Molecule_Category
from ctgomartini.core.Molecule import Bonded_Categories, Molecule_Categories

# ============================================
# 示例 1: 创建一个简单的分子
# ============================================
print("=" * 60)
print("示例 1: 创建分子并添加原子")
print("=" * 60)

# 创建一个名为 "Protein" 的分子
protein = Molecule(name="Protein")
print(f"创建分子: {protein.name}")

# 模拟从拓扑文件读取 atoms 节的数据
# 格式: atomnr, atomtype, resnr, resname, atomname, cgnr, charge, mass
atom_lines = [
    "1 Q5 1 ALA BB 1 0.0 72.0",
    "2 SC2 1 ALA SC1 2 0.0 54.0", 
    "3 P2 2 GLY BB 3 0.0 72.0",
]

for line in atom_lines:
    fields = line.split()
    protein.readLine(line, "atoms")  # 分类存储到 atoms 节

print(f"已添加 {len(protein._topology['atoms'])} 个原子")
for atom in protein._topology['atoms']:
    print(f"  Atom {atom[0]}: {atom[4]} ({atom[1]}), Res {atom[2]} {atom[3]}")


# ============================================
# 示例 2: 添加键合相互作用
# ============================================
print("\n" + "=" * 60)
print("示例 2: 添加键合相互作用")
print("=" * 60)

# 添加 bonds (格式: atom1, atom2, functype, length, force_constant)
bond_lines = [
    "1 2 1 0.35 5000.0",  # BB-SC1 键
    "2 3 1 0.38 5000.0",  # SC1-BB 键
]

for line in bond_lines:
    protein.readLine(line, "bonds")

print(f"已添加 {len(protein._topology['bonds'])} 个键")
for bond in protein._topology['bonds']:
    print(f"  Bond: atom {bond[0]} - atom {bond[1]}, length={bond[3]} nm, k={bond[4]} kJ/mol/nm²")


# ============================================
# 示例 3: 查看 Bonded_Categories 注册表
# ============================================
print("\n" + "=" * 60)
print("示例 3: Bonded_Categories 注册表")
print("=" * 60)

print("已注册的键合相互作用类别:")
for category_name in Bonded_Categories:
    cat = Bonded_Categories[category_name]
    print(f"  - {category_name}: {cat.description}")


# ============================================
# 示例 4: 创建自定义分子类别
# ============================================
print("\n" + "=" * 60)
print("示例 4: 创建自定义分子类别")
print("=" * 60)

# 创建一个自定义类别（例如用于存储特殊参数）
custom_category = Molecule_Category(
    name="custom_params",
    description="Custom parameters for analysis",
    category="custom_params",
    contents=[
        ["param1", "10.0"],
        ["param2", "20.0"],
    ]
)

print(f"创建类别: {custom_category.name}")
print(f"描述: {custom_category.description}")
print(f"内容: {custom_category.contents}")


# ============================================
# 示例 5: 模拟从 Topology 解析的完整流程
# ============================================
print("\n" + "=" * 60)
print("示例 5: 完整解析流程模拟")
print("=" * 60)

# 创建一个新的多势阱蛋白分子
mbp_protein = Molecule(name="GlnBP_MBP")

# 1. 添加 moleculetype
mbp_protein.readLine("GlnBP_MBP 1", "moleculetype")

# 2. 添加 multiple_basin 参数
mbp_protein.readLine("True exp 2 1/1000 -50 0", "multiple_basin")

# 3. 添加多态角度（不同态有不同的角度值）
# 格式: atom1 atom2 atom3 n_states state_id functype angle k
multi_angle_lines = [
    "10 11 12 2 1 10 127.0 20.0",  # State 1: Open 构象
    "10 11 12 2 2 10 100.0 20.0",  # State 2: Closed 构象
]

for line in multi_angle_lines:
    mbp_protein.readLine(line, "multi_angles")

# 查看存储的数据结构
print(f"\n分子名称: {mbp_protein._topology['moleculetype'][0][0]}")
print(f"MBP 参数: {mbp_protein._topology['multiple_basin']}")
print(f"Multi-angles 条目数: {len(mbp_protein._topology['multi_angles'])}")

# 显示 multi_angles 详情
print("\n多态角度定义:")
for angle in mbp_protein._topology['multi_angles']:
    state_id = angle[4]
    angle_val = angle[6]
    print(f"  State {state_id}: {angle_val}°")


# ============================================
# 示例 6: Molecule_Category 的补充参数功能
# ============================================
print("\n" + "=" * 60)
print("示例 6: Atoms 类别自动补充参数")
print("=" * 60)

# 创建 Atoms 类别（带参数补充功能）
atoms_cat = Molecule_Categories['atoms']([
    ["1", "Q5", "1", "ALA", "BB", "1", None, None],  # 缺少 charge 和 mass
])

# 模拟力场参数 (格式: [name, bonded_type, atomic_num, mass, charge, particle_type, V, W])
ff_params = {
    'atomtypes': [
        ['Q5', None, None, '72.0', '0.0', 'A', '0.0', '0.0']
    ]
}

print("补充前:")
print(f"  Charge: {atoms_cat.contents[0][6]}")
print(f"  Mass: {atoms_cat.contents[0][7]}")

# 调用 complement 补充缺失参数
atoms_cat.complement(ff_params)

print("\n补充后:")
print(f"  Charge: {atoms_cat.contents[0][6]}")
print(f"  Mass: {atoms_cat.contents[0][7]}")


print("\n" + "=" * 60)
print("示例完成!")
print("=" * 60)
