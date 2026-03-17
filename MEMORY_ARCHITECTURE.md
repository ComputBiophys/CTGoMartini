# CTGoMartini 架构设计文档

> 生成时间: 2026-03-17  
> 用途: REMD功能整合后的架构方案

---

## 一、目录架构

```
ctgomartini/
├── __init__.py
├── _version.py
│
├── ctgomartinize.py        # 拓扑生成主入口
├── run_ctgomartini.py      # 统一运行入口（自动识别MD/REMD）
│
├── api/                    # 用户API
│   ├── __init__.py
│   ├── MartiniTopology.py  # MartiniTopFile类
│   ├── MBMoleculeTop.py    # GenMBPTop多盆地拓扑生成
│   └── SBMoleculeTop.py    # GenSBPTop单盆地拓扑生成
│
├── core/                   # 核心（数据 + 模拟）
│   ├── __init__.py
│   ├── Topology.py         # 纯拓扑解析（GROMACS格式）
│   ├── Molecule.py         # 分子表示
│   ├── ForceField.py       # 力场参数
│   ├── bonded/             # 键合相互作用类型
│   │   ├── __init__.py
│   │   ├── bonds.py
│   │   ├── angles.py
│   │   ├── dihedrals.py
│   │   ├── virtual_sites.py
│   │   └── multi_state.py  # 多盆地相互作用
│   ├── Nonbonded_interaction.py
│   ├── CombineMols.py
│   ├── vsites.py
│   │
│   └── simulation/         # 【整合进core】模拟运行核心
│       ├── __init__.py     # 导出 Runner, MDRunner, REMDRunner
│       ├── base.py         # BaseRunner: 平台加载、结构读取、约束、压浴
│       ├── md.py           # MDRunner(BaseRunner): 标准MD运行
│       ├── remd.py         # REMDRunner(BaseRunner): REMD + 重启功能
│       └── config.py       # 输入配置解析（识别REMD参数）
│
├── utils/                  # 通用工具（瘦身）
│   ├── __init__.py
│   ├── ReadInp.py          # 逐步迁移到core/simulation/config.py
│   ├── WriteItp.py         # ITP文件写入
│   ├── constraints_to_bonds.py
│   └── ...
│
├── data/                   # 【纯数据】力场、模板、示例
│   ├── ForceFields/
│   ├── Membrane/
│   ├── Soluble/
│   └── OVrCSU.py           # 可选：辅助脚本
│
└── analysis/               # 分析工具
    ├── __init__.py
    └── ...
```

---

## 二、pyproject.toml 配置

```toml
[project.scripts]
ctgomartinize = "ctgomartini.ctgomartinize:main"
run_ctgomartini = "ctgomartini.run_ctgomartini:main"
```

---

## 三、使用方式

### 1. 拓扑生成

```bash
ctgomartinize -s state1.pdb state2.pdb -m map1.map map2.map -mol Open Closed -mbmol GlnBP
```

### 2. 模拟运行（自动识别）

```bash
# 普通MD（inp中无exc_freq参数）
run_ctgomartini -i md.inp

# REMD（inp中有exc_freq参数，自动识别）
run_ctgomartini -i md_REMD.inp

# REMD重启（检测已有output.nc或指定--restart）
run_ctgomartini -i md_REMD.inp --restart
```

---

## 四、inp文件区分

### 普通MD (`md.inp`)

```ini
nstep       = 50000000
dt          = 0.020
temp        = 310
pcouple     = yes
# ... 无exc_freq参数
```

### REMD (`md_REMD.inp`)

```ini
nstep       = 50000000
dt          = 0.020
temp        = 310
pcouple     = yes
exc_freq    = 1000        # ← REMD特有参数，自动识别为REMD模式
```

---

## 五、运行流程

```
run_ctgomartini.py -i xxx.inp
    │
    ├── 1. 解析命令行参数 (-i, --restart)
    │
    ├── 2. 读取并解析inp文件
    │   └── core/simulation/config.py
    │
    ├── 3. 检测exc_freq参数？
    │       ├── 不存在 → 创建MDRunner → run()
    │       └── 存在 → 创建REMDRunner
    │               ├── 检测--restart或output.nc存在？
    │               │       ├── 否 → run()      # 新REMD
    │               │       └── 是 → extend()   # 重启续跑
    │               └── 执行
    │
    └── 4. 执行对应运行器
```

---

## 六、各模块职责

| 模块 | 职责 |
|------|------|
| `ctgomartinize.py` | **CLI入口**：解析参数 → 调用`api/`生成拓扑 |
| `run_ctgomartini.py` | **统一CLI入口**：解析inp → 自动选择运行器 |
| `api/` | **用户API**：MartiniTopFile, GenMBPTop, GenSBPTop |
| `core/` | **数据+模拟**：根级数据结构 + `simulation/`运行控制 |
| `core/simulation/base.py` | **BaseRunner**：平台加载、结构读取、约束、压浴等公共逻辑 |
| `core/simulation/md.py` | **MDRunner**：标准MD运行 |
| `core/simulation/remd.py` | **REMDRunner**：REMD运行 + 重启续跑 |
| `core/simulation/config.py` | **配置解析**：识别REMD参数（exc_freq）、重启状态 |
| `utils/` | **通用工具**：I/O、格式转换、辅助函数 |
| `data/` | **纯数据**：力场文件、模板、示例（无Python业务代码） |
| `analysis/` | **分析工具**：后处理分析 |

---

## 七、伪代码示意

### run_ctgomartini.py

```python
#!/usr/bin/env python3
"""统一运行入口：自动识别MD/REMD模式"""

import argparse
from ctgomartini.core.simulation import MDRunner, REMDRunner
from ctgomartini.core.simulation.config import load_config

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input parameter file")
    parser.add_argument("--restart", action="store_true", help="Restart from checkpoint")
    args = parser.parse_args()
    
    # 解析配置文件
    config = load_config(args.input)
    
    # 根据配置自动选择运行模式
    if config.exc_freq is not None:
        # REMD模式
        runner = REMDRunner(config)
        if args.restart or output_nc_exists(config.output_data):
            runner.extend()   # 重启续跑
        else:
            runner.run()      # 新REMD模拟
    else:
        # 普通MD模式
        runner = MDRunner(config)
        runner.run()

if __name__ == "__main__":
    main()
```

### core/simulation/config.py

```python
"""输入配置解析"""

class SimulationConfig:
    """模拟配置类"""
    def __init__(self):
        # 通用参数
        self.nstep: int = 0
        self.dt: float = 0.020
        self.temp: float = 310.0
        
        # REMD特有参数（None表示非REMD）
        self.exc_freq: int | None = None
        
        # 文件路径
        self.input: str = "input.gro"
        self.topol: str = "system.top"
        self.output_data: str = "output.nc"
        
        # ... 其他参数

def load_config(inpfile: str) -> SimulationConfig:
    """从inp文件加载配置"""
    config = SimulationConfig()
    # 解析GROMACS风格参数文件
    # 检测exc_freq参数是否存在
    return config
```

### core/simulation/base.py

```python
"""运行器基类"""

import openmm as mm
from ctgomartini.api import MartiniTopFile

class BaseRunner:
    """模拟运行器基类"""
    
    def __init__(self, config: SimulationConfig):
        self.config = config
        self.platform = None
        self.simulation = None
    
    def load_platform(self):
        """加载OpenMM平台（CPU/CUDA/OpenCL）"""
        pass
    
    def load_structure(self):
        """加载结构文件（GRO/PDB）"""
        pass
    
    def create_system(self):
        """创建OpenMM System"""
        pass
    
    def add_restraints(self):
        """添加位置约束"""
        pass
    
    def add_barostat(self):
        """添加压浴"""
        pass
```

### core/simulation/md.py

```python
"""标准MD运行器"""

from .base import BaseRunner

class MDRunner(BaseRunner):
    """标准分子动力学运行器"""
    
    def run(self):
        """执行MD模拟"""
        # 1. 设置平台
        # 2. 加载结构
        # 3. 创建系统
        # 4. 添加约束/压浴
        # 5. 创建积分器
        # 6. 运行模拟循环
        pass
```

### core/simulation/remd.py

```python
"""REMD运行器"""

import openmmtools
from openmmtools.multistate import ReplicaExchangeSampler, MultiStateReporter
from .base import BaseRunner

class REMDRunner(BaseRunner):
    """副本交换分子动力学运行器"""
    
    def __init__(self, config: SimulationConfig):
        super().__init__(config)
        self.replica_params = self._build_replica_params()
    
    def _build_replica_params(self):
        """构建副本参数（beta, C1, C2数组）"""
        # 从config解析或生成线性间距的C1参数
        pass
    
    def run(self):
        """执行新REMD模拟"""
        # 1. 为每个副本创建不同参数的system
        # 2. 创建thermodynamic_states
        # 3. 创建ReplicaExchangeSampler
        # 4. 运行模拟
        pass
    
    def extend(self, n_iterations: int | None = None):
        """从检查点续跑"""
        # 1. 从output.nc加载检查点
        # 2. 计算需要补充的步数
        # 3. 继续运行
        pass
```

---

## 八、关键设计决策

### 1. 自动模式识别
- 通过检测 `exc_freq` 参数自动区分MD/REMD
- 无需用户手动指定 `--mode`

### 2. simulation整合到core
- `core/` 包含：数据结构 + 运行控制
- 保持顶层模块数量精简（api, core, utils, data, analysis）

### 3. data目录纯数据化
- 原 `data/run_REMD.py` 等代码移至 `core/simulation/`
- `data/` 仅保留力场、模板、示例文件

### 4. 统一入口
- 一个 `run_ctgomartini` 命令处理所有模拟类型
- 重启功能通过 `--restart` 标志或自动检测实现

---

## 九、依赖关系

```
ctgomartinize.py
    └── api/ (MartiniTopFile, GenMBPTop, GenSBPTop)
    └── utils/ (write_itp, read_inputs)

run_ctgomartini.py
    └── core/simulation/ (MDRunner, REMDRunner)
    └── core/simulation/config.py (load_config)
    
core/simulation/md.py
    └── core/simulation/base.py (BaseRunner)
    └── api/ (MartiniTopFile)
    
core/simulation/remd.py
    └── core/simulation/base.py (BaseRunner)
    └── api/ (MartiniTopFile)
    └── openmmtools (ReplicaExchangeSampler)
```

---

## 十、待完成任务清单

- [ ] 创建 `core/simulation/` 目录
- [ ] 将 `data/Topology.py` 内容迁移至 `core/simulation/base.py`
- [ ] 重构 `data/run_ctgomartini.py` → `run_ctgomartini.py`（包根）
- [ ] 重构 `data/run_REMD.py` + `data/run_REMD_restart.py` → `core/simulation/remd.py`
- [ ] 创建 `core/simulation/config.py` 统一配置解析
- [ ] 更新 `pyproject.toml` 的 `[project.scripts]`
- [ ] 清理 `data/` 目录，移除非数据文件
- [ ] 更新 `utils/ReadInp.py`，逐步迁移功能到 `core/simulation/config.py`
