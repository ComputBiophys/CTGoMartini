# CTGoMartini 重构状态记录

## 重构目标

将 CTGoMartini 代码库重构为模块化、现代化的架构：
- 统一拓扑模块：`api/` + `core/` → `topology/`
- 新增仿真模块：`simulation/`
- 新增 CLI 模块：`cli/`
- 清理数据目录：`data/` 只保留数据和模板
- 命名规范：camelCase/PascalCase → snake_case

---

## 重构完成状态 (2026-03-18)

### ✅ 1. 拓扑模块重组 (完成)

| 原路径 | 新路径 | 状态 |
|--------|--------|------|
| `ctgomartini/api/MartiniTopology.py` | `ctgomartini/topology/builder.py` | ✅ MartiniTopFile |
| `ctgomartini/api/MBMoleculeTop.py` | `ctgomartini/topology/generator/multi_basin.py` | ✅ create_mb_topology |
| `ctgomartini/api/SBMoleculeTop.py` | `ctgomartini/topology/generator/single_basin.py` | ✅ create_sb_topology |
| `ctgomartini/core/Topology.py` | `ctgomartini/topology/parser.py` | ✅ TopologyParser |
| `ctgomartini/core/Molecule.py` | `ctgomartini/topology/models.py` | ✅ Molecule, ForceField |
| `ctgomartini/core/ForceField.py` | `ctgomartini/topology/models.py` | ✅ 合并到 models.py |
| `ctgomartini/core/CombineMols.py` | `ctgomartini/topology/generator/combiner.py` | ✅ combine_* 函数 |
| `ctgomartini/core/bonded/` | `ctgomartini/topology/interactions/` | ✅ 所有相互作用类 |
| `ctgomartini/core/Nonbonded_interaction.py` | `ctgomartini/topology/interactions/nonbonded.py` | ✅ |
| `ctgomartini/core/vsites.py` | `ctgomartini/topology/interactions/vsites.py` | ✅ |

### ✅ 2. 仿真模块创建 (完成)

| 新文件 | 功能 | 状态 |
|--------|------|------|
| `ctgomartini/simulation/__init__.py` | 模块导出 | ✅ |
| `ctgomartini/simulation/config.py` | SimulationConfig, load_config | ✅ |
| `ctgomartini/simulation/base.py` | SimulationRunner 基类 | ✅ |
| `ctgomartini/simulation/md.py` | MDRunner 标准 MD | ✅ |
| `ctgomartini/simulation/remd.py` | REMDRunner 副本交换 | ✅ |

### ✅ 3. CLI 模块创建 (完成)

| 新文件 | 功能 | 状态 |
|--------|------|------|
| `ctgomartini/cli/__init__.py` | 模块导出 | ✅ |
| `ctgomartini/cli/run_ctgomartini.py` | 统一仿真运行器 | ✅ |
| `ctgomartini/cli/ctgomartinize.py` | 拓扑生成器 | ✅ |

**入口点** (pyproject.toml):
- `run_ctgomartini` → `ctgomartini.cli.run_ctgomartini:main`
- `ctgomartinize` → `ctgomartini.cli.ctgomartinize:main`

### ✅ 4. 工具模块重组 (完成)

| 原文件 | 新路径 | 状态 |
|--------|--------|------|
| `utils/WriteItp.py` | `utils/write_itp.py` | ✅ 蛇形命名 |
| `utils/ReadInp.py` | **已删除** | ✅ 功能在 simulation/config.py |
| `utils/ConvertLongShortElasticBonds.py` | `utils/convert_long_short_elastic_bonds.py` | ✅ |
| `utils/Create_goVirt_for_multimer.py` | **已删除** | ✅ 被 martinize2 -go 替代 |
| `utils/OVrCSU.py` | `utils/contact_map.py` | ✅ 重命名 |
| 新增 | `utils/constraints_to_bonds.py` | ✅ 约束转键 |
| 新增 | `utils/pdb_validation.py` | ✅ PDB 验证 |

### ✅ 5. 数据目录清理 (完成)

**删除的文件**:
- `data/run_ctgomartini.py` → 被 cli/ 替代
- `data/run_REMD.py` → 被 simulation/remd.py 替代
- `data/run_REMD_restart.py` → 被 simulation/remd.py 替代
- `data/Topology.py` → 功能在 simulation/base.py
- `data/Constraint2Bonds.py` → 功能在 utils/constraints_to_bonds.py
- `data/ctgomartinize.ipynb` → 删除
- `data/REMD/*.ipynb` (5个) → 转化为 analysis/remd_*.py

**保留的目录**:
- `data/ForceFields/` - 力场文件
- `data/Membrane/` - 膜蛋白模板
- `data/Soluble/` - 可溶性蛋白模板

### ✅ 6. 分析模块重组 (完成)

| 原文件 | 新路径 | 状态 |
|--------|--------|------|
| `analysis/MBAnalysis_v2_argparse.py` | **已删除** | ✅ 旧版本 |
| `analysis/MBAnalysis_v3_argparse.py` | **已删除** | ✅ 旧版本 |
| `analysis/MBAnalysis_v4_argparse.py` | `analysis/drms_analysis.py` | ✅ 重命名 |
| `analysis/QValueAnalysis.ipynb` | **已删除** | ✅ |
| `data/REMD/ExchangeRatio.py` | `analysis/remd_exchange_ratio.py` | ✅ |
| 新增 | `analysis/remd_drms_analysis.py` | ✅ 从 notebook 转换 |
| 新增 | `analysis/remd_replica_state.py` | ✅ 从 notebook 转换 |
| 新增 | `analysis/remd_free_energy.py` | ✅ 从 notebook 转换 |
| 新增 | `analysis/remd_mbar_analysis.py` | ✅ 从 notebook 转换 |

---

## 命名规范转换

| 旧名称 (PascalCase) | 新名称 (snake_case) | 位置 |
|---------------------|---------------------|------|
| `Topology` | `TopologyParser` | topology.parser |
| `GenMBPTop` | `create_mb_topology` | topology.generator |
| `GenSBPTop` | `create_sb_topology` | topology.generator |
| `read_inputs` | `load_config` | simulation.config |
| `WriteItp` | `write_itp` | utils.write_itp |
| `ConvertLongShortElasticBonds` | `convert_long_short_elastic_bonds` | utils |

---

## 当前测试状态 (2026-03-18)

### 通过测试 (137/137, 100%)

| 测试类别 | 数量 | 状态 |
|----------|------|------|
| `tests/unit/` | 84 | ✅ 全部通过 |
| `tests/integration/` | 44 | ✅ 全部通过 |
| `tests/e2e/` | 9 | ✅ 全部通过 |
| **总计** | **137** | **✅ 100%** |

### 关键测试覆盖

- ✅ Classic Martini 测试 (30个)
- ✅ Contacts 测试 (8个)
- ✅ Multiple Basin Potential 测试 (4个)
- ✅ Energy Item Comparison 测试 (2个)
- ✅ ITP 生成功能测试
- ✅ 约束转键测试
- ✅ PDB 验证测试
- ✅ Vermouth 0.15+ 兼容性测试 (5个E2E)

---

## 文档更新状态

| 文档 | 状态 | 说明 |
|------|------|------|
| `README.md` | ✅ 已更新 | 新模块结构、使用示例 |
| `AGENTS.md` | ✅ 已更新 | 项目结构、导入示例 |
| `MIGRATION_GUIDE.md` | ✅ 已更新 | 1.0.0 迁移指南 |
| `REFACTORING_STATUS.md` | ✅ 已更新 | 本文件 |

---

## 导入示例

### 拓扑模块

```python
from ctgomartini.topology import (
    TopologyParser,          # GROMACS 拓扑解析
    MartiniTopFile,          # OpenMM 系统构建
    Molecule,                # 分子模型
    ForceField,              # 力场参数
    create_sb_topology,      # 单盆地拓扑生成
    create_mb_topology,      # 多盆地拓扑生成
    combine_atoms,           # 分子组合函数
    combine_bonds_constraints,
    combine_angles,
    combine_dihedrals,
    combine_contacts,
)

# 相互作用类
from ctgomartini.topology.interactions import (
    HarmonicBonds,
    ContactsLJ, Contacts6_12, Contacts10_12,
    HarmonicAngles, G96Angles, RestrictedAngles,
    PeriodicDihedrals,
    EXPInteraction, HAMInteraction,
)
```

### 仿真模块

```python
from ctgomartini.simulation import (
    SimulationConfig,    # 配置类
    load_config,         # 加载 .inp 文件
    SimulationRunner,    # 基类
    MDRunner,            # 标准 MD
    REMDRunner,          # 副本交换 MD
)
```

### 工具模块

```python
from ctgomartini.utils import (
    write_itp,
    convert_map_format,
    convert_long_short_elastic_bonds,
    convert_constraints_to_bonds,
    validate_pdb_compatibility,
)
```

### 分析模块

```python
from ctgomartini.analysis import NativeContanceAnalysis
from ctgomartini.analysis.drms_analysis import DRMSAnalyzer

# REMD 分析 (需要 openmmtools)
from ctgomartini.analysis.remd_drms_analysis import plot_drms_trajectory
from ctgomartini.analysis.remd_replica_state import plot_replica_states
```

---

## 环境配置建议

### 当前推荐环境

```bash
# Python 环境
python = 3.12+

# 核心依赖
openmm >= 8.1
mdanalysis >= 2.4
vermouth >= 0.15.0
numpy, pandas, matplotlib

# 可选依赖
openmm-plumed >= 2.0.0  # 增强采样
pymbar                   # REMD 分析
openmmtools              # REMD 仿真
```

---

## 重构完成总结

### 主要成就

1. **✅ 模块结构清晰化** - 从扁平结构到分层架构
2. **✅ 命名规范化** - 统一使用 snake_case
3. **✅ 功能分离** - 拓扑、仿真、CLI、分析各自独立
4. **✅ 向后兼容** - 提供清晰的迁移指南
5. **✅ 测试覆盖** - 137个测试全部通过
6. **✅ 文档同步** - 所有文档已更新

### 目录结构对比

```
重构前 (0.6.x):
ctgomartini/
├── api/          # 高层 API (3个文件)
├── core/         # 核心功能 (10+个文件)
├── util/         # 工具函数 (旧命名)
├── data/         # 数据和脚本 (混乱)
└── analysis/     # 分析脚本 (不完整)

重构后 (1.0.0):
ctgomartini/
├── topology/     # 统一拓扑模块 (清晰分层)
│   ├── parser.py
│   ├── models.py
│   ├── builder.py
│   ├── generator/
│   └── interactions/
├── simulation/   # 仿真执行模块
│   ├── config.py
│   ├── base.py
│   ├── md.py
│   └── remd.py
├── cli/          # 命令行接口
│   ├── run_ctgomartini.py
│   └── ctgomartinize.py
├── utils/        # 工具函数 (蛇形命名)
├── analysis/     # 完整分析套件
│   ├── QValueAnalysis.py
│   ├── drms_analysis.py
│   └── remd_*.py
└── data/         # 仅数据和模板
    ├── ForceFields/
    ├── Membrane/
    └── Soluble/
```

---

**记录时间**: 2026-03-18  
**重构版本**: 0.6.x → 1.0.0  
**测试通过率**: 100% (137/137) ✅  
**状态**: **重构完成**
