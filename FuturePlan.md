================================================================================
           新架构：5大功能模块 + 规范命名
================================================================================

ctgomartini/
│
├── 📁 topology/                    # 【功能1: 拓扑构建】
│   ├── __init__.py
│   │   └── 暴露: TopologyParser, TopologyBuilder, 
│   │           MartiniTopFile, create_mb_topology, create_sb_topology
│   │
│   ├── parser.py                   # 拓扑解析
│   │   └── class TopologyParser    # 原 Topology
│   │   └── class TopologyData      # 原 Molecule + ForceField 整合
│   │
│   ├── builder.py                  # 拓扑构建器
│   │   └── class TopologyBuilder   # 原 MartiniTopFile
│   │   └── def create_mb_topology  # 原 GenMBPTop
│   │   └── def create_sb_topology  # 原 GenSBPTop
│   │
│   ├── combiner.py                 # 分子组合
│   │   └── def combine_molecules   # 原 CombineTwoMolecules
│   │   └── class CombinedTopology  # 原 CombinedMolecule
│   │
│   ├── virtual_sites.py            # 虚拟位点
│   │   └── class VirtualSiteBuilder # 原 VirtualSites
│   │
│   └── interactions/               # 相互作用定义
│       ├── __init__.py
│       ├── base.py                 # 基础类
│       │   └── class Interaction
│       │   └── class InteractionError
│       ├── registry.py             # 注册表
│       │   └── def register_interaction
│       │   └── def get_interaction_handler
│       ├── bonds.py                # 键
│       │   └── class HarmonicBond
│       │   └── class ConstraintBond
│       ├── angles.py               # 角
│       │   └── class G96Angle
│       │   └── class RestrictedAngle
│       ├── dihedrals.py            # 二面角
│       │   └── class PeriodicDihedral
│       │   └── class RyckaertBellemansDihedral
│       ├── contacts.py             # 接触势
│       │   └── class LJContact
│       │   └── class AdjustableLJContact
│       │   └── class LJ1012Contact
│       ├── nonbonded.py            # 非键作用
│       │   └── class NonbondedForce
│       │   └── class ElectrostaticCorrection
│       ├── multi_state.py          # 多状态交互
│       │   └── class MultiStateInteraction
│       │   └── class MultiAllBonds
│       └── mixing.py               # 混合方案（已修改）
│           └── class EXPMixing      # 原 EXPInteraction
│           └── class HAMMixing      # 原 HAMInteraction
│
├── 📁 simulation/                  # 【功能2: 分子动力学模拟】
│   ├── __init__.py
│   │   └── 暴露: SimulationRunner, MDRunner, REMDRunner, 
│   │           load_config, SimulationConfig
│   │
│   ├── config.py                   # 配置解析
│   │   └── class SimulationConfig  # 原 _OpenMMReadInputs
│   │   └── def load_config         # 原 read_inputs
│   │
│   ├── base.py                     # 模拟基类
│   │   └── class SimulationRunner  # 原 BaseRunner
│   │   └── def write_checkpoint
│   │   └── def add_restraints
│   │   └── def generate_restraints # 原 generate_restraints
│   │
│   ├── md.py                       # 标准MD
│   │   └── class MDRunner          # 保持
│   │   └── def run_md_simulation   # 辅助函数
│   │
│   ├── remd.py                     # 副本交换MD
│   │   └── class REMDRunner        # 保持
│   │   └── def run_remd_simulation
│   │
│   └── reporters.py                # 输出报告（新）
│       └── class CheckpointReporter
│       └── class EnergyReporter
│
├── 📁 cli/                         # 【功能3: 命令行接口】
│   ├── __init__.py
│   ├── main.py                     # 主入口
│   │   └── def main                # 统一入口
│   ├── martinize.py                # 拓扑生成命令
│   │   └── def martinize_cmd       # 原 ctgomartinize.py
│   ├── run.py                      # 运行命令
│   │   └── def run_cmd             # 原 run_ctgomartini.py
│   └── analyze.py                  # 分析命令（新）
│       └── def analyze_cmd
│
├── 📁 analysis/                    # 【功能4: 分析工具】
│   ├── __init__.py
│   │   └── 暴露: QValueAnalyzer, MultiBasinAnalyzer
│   │
│   ├── qvalue.py                   # Q值分析
│   │   └── class QValueAnalyzer    # 原 QValueAnalysis 中的功能
│   │   └── def compute_q_values    # 原 computeQValue
│   │   └── def analyze_trajectory
│   │
│   ├── multibasin.py               # 多盆地分析
│   │   └── class MultiBasinAnalyzer # 整合 MBAnalysis_v*
│   │   └── def analyze_state_distribution
│   │   └── def compute_free_energy
│   │
│   └── trajectory.py               # 轨迹工具
│       └── def load_trajectory     # 统一加载
│       └── def align_trajectory
│
├── 📁 utils/                       # 【功能5: 通用工具】
│   ├── __init__.py
│   │
│   ├── io/                         # I/O相关
│   │   ├── __init__.py
│   │   ├── itp.py                  # ITP文件
│   │   │   └── def write_itp       # 保持
│   │   │   └── def read_itp
│   │   ├── pdb.py                  # PDB/GRO处理
│   │   │   └── def load_structure  # 原 load_structure
│   │   │   └── def write_pdb
│   │   │   └── def write_gro
│   │   └── formats.py              # 格式转换
│   │       └── def convert_map_format
│   │       └── def convert_elastic_bonds
│   │
│   ├── contacts/                   # 接触映射工具
│   │   ├── __init__.py
│   │   ├── rcsu.py                 # rCSU相关
│   │   │   └── def generate_contact_map
│   │   └── ov.py                   # OV方法
│   │       └── def ov_contact_map
│   │
│   ├── multimer.py                 # 多聚体工具
│   │   └── def create_virtual_sites_for_multimer
│   │   └── def setup_go_multimer
│   │
│   └── validation.py               # 验证工具
│       └── def validate_pdb
│       └── def check_topology
│
├── 📁 data/                        # 【数据文件】（简化）
│   ├── forcefields/                # 原 ForceFields/
│   │   └── martini300/
│   ├── templates/                  # 原 Membrane/, Soluble/
│   │   ├── membrane/
│   │   └── soluble/
│   └── notebooks/                  # 原 REMD/ 的 ipynb
│       └── README.md
│
├── __init__.py                     # 包初始化
│   └── 暴露主要接口
├── __main__.py                     # python -m ctgomartini
│   └── from cli.main import main; main()
└── _version.py                     # 版本信息

================================================================================
                        命名对照表
================================================================================

【类名 PascalCase】
原名称                          新名称                          位置
────────────────────────────────────────────────────────────────────────────
Topology                        TopologyParser                  topology/parser.py
Molecule + ForceField           TopologyData                    topology/parser.py
MartiniTopFile                  TopologyBuilder                 topology/builder.py
GenMBPTop (函数)                create_mb_topology              topology/builder.py
GenSBPTop (函数)                create_sb_topology              topology/builder.py
CombineTwoMolecules (函数)      combine_molecules               topology/combiner.py
CombinedMolecule                CombinedTopology                topology/combiner.py
VirtualSites                    VirtualSiteBuilder              topology/virtual_sites.py
OpenMMReadInputs                SimulationConfig                simulation/config.py
BaseRunner                      SimulationRunner                simulation/base.py
MDRunner                        MDRunner                        simulation/md.py (不变)
REMDRunner                      REMDRunner                      simulation/remd.py (不变)
EXPInteraction                  EXPMixing                       topology/interactions/mixing.py
HAMInteraction                  HAMMixing                       topology/interactions/mixing.py
Nonbonded_interaction           NonbondedForce                  topology/interactions/nonbonded.py
Interaction (基础)              Interaction                     topology/interactions/base.py

【函数名 snake_case】
原名称                          新名称                          位置
────────────────────────────────────────────────────────────────────────────
read_inputs                     load_config                     simulation/config.py
write_itp                       write_itp                       utils/io/itp.py (不变)
load_structure                  load_structure                  utils/io/pdb.py (不变)
convert_map_format              convert_map_format              utils/io/formats.py (不变)
generate_restraints             generate_restraints             simulation/base.py (不变)

================================================================================
                        用户代码示例
================================================================================

# 1. 拓扑构建
from ctgomartini.topology import TopologyBuilder, create_mb_topology

builder = TopologyBuilder("system.top")
system = builder.build()

# 或使用函数
system = create_mb_topology(
    state_files=["closed.itp", "open.itp"],
    method="EXP",
    beta=0.00333,
    offsets=[-300.0, 0.0]
)

# 2. 模拟
from ctgomartini.simulation import MDRunner, load_config

config = load_config("md.inp")
runner = MDRunner(config)
runner.run()

# 3. 命令行
# ctgomartini martinize -s state1.pdb state2.pdb ...
# ctgomartini run -i md.inp
# ctgomartini analyze -t traj.xtc --qvalue

# 4. 分析
from ctgomartini.analysis import QValueAnalyzer

analyzer = QValueAnalyzer(topology="system.top", trajectory="traj.xtc")
results = analyzer.compute_q_values()