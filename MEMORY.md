# CTGoMartini 项目关键记忆

> 本文件记录项目重构的关键决策和状态，用于快速恢复上下文。
> 最后更新: 2026-03-13

## 一、项目架构（当前）

```
CTGoMartini/
├── ctgomartini/          # 主包（重构后）
│   ├── api/              # MartiniTopFile, GenMBPTop, GenSBPTop
│   ├── core/             # 原 func/ - Topology, Molecule, ForceField
│   ├── utils/            # 原 util/ - write_itp, read_inputs
│   ├── data/             # 可执行脚本
│   └── analysis/         # QValueAnalysis
├── tests/                # 测试套件（已重构）
│   ├── unit/             # 原 func/ 测试
│   ├── integration/      # 原 api/ 测试
│   ├── e2e/              # 原 run/ 测试
│   └── fixtures/         # 原 data/ 测试数据
└── Tutorial/             # GlnBP, TRAAK 示例
```

## 二、命名规范变更

| 旧命名 | 新命名 | 类型 |
|--------|--------|------|
| `func/` | `core/` | 目录 |
| `util/` | `utils/` | 目录 |
| `WriteItp` | `write_itp` | 函数 |
| `ConvertLongShortElasticBonds` | `convert_long_short_elastic_bonds` | 函数 |
| `Create_goVirt_for_multimer` | `create_go_virt_for_multimer` | 函数 |
| `gen_restraints` | `generate_restraints` | 函数 |
| `restraints` | `add_restraints` | 函数 |

## 三、关键 Bug 修复记录

### 1. Local_BondedInteraction_dict 导入错误
- **问题**: 从 Nonbonded_interaction 导入了空字典
- **修复**: 改为从 Bonded_interaction 导入
- **文件**: `ctgomartini/core/__init__.py`

### 2. MBP 能量表达式为空
- **问题**: `energy1 = ; energy2 = ;` 因为循环遍历了错误的类别
- **修复**: 修改循环为固定的 `['multi_angles', 'multi_dihedrals', 'multi_contacts']`
- **文件**: `ctgomartini/api/MartiniTopology.py:264`

### 3. 竞态条件（多进程写入）
- **问题**: `test_MBGoMartini.py` 多进程同时写入同一 ITP 文件
- **修复**: 使用 `WorkingDirectoryContext` 上下文管理器，并将 `SetMBPParameter` 移到多进程外
- **文件**: `tests/integration/test_MBGoMartini.py`

### 4. DSSP 版本兼容性
- **问题**: vermouth 0.9.6 只支持 DSSP 2.2.1/3.0.0
- **修复**: 测试优先使用 `mkdssp`，主脚本移除 `-dssp` 参数依赖

## 四、测试架构

| 类别 | 数量 | 位置 | 说明 |
|------|------|------|------|
| unit | 2 | `tests/unit/` | ITP 写入、键转换 |
| integration | 42 | `tests/integration/` | OpenMM/GROMACS 对比、MBP |
| e2e | 6 | `tests/e2e/` | 完整工作流 |

**测试状态**: 50/50 通过 (100%)

## 五、待办事项（按优先级）

### 🔴 高优先级

1. **Vermouth 0.15.0 升级**
   - Contact Maps: 保留 OV+rCSU，添加 `--use-martinize-go` 选项
   - Secondary Structure: 使用 MDTraj，移除 DSSP
   - CLI 参数: 显式指定 -scfix on, -ef 700, -cys auto
   - 时间: 2 周

### 🟡 中优先级

2. **文档完善**
   - API 参考文档 (Sphinx/MkDocs)
   - 教程更新

3. **CI/CD**
   - GitHub Actions 自动化测试

### 🟢 低优先级

4. **功能增强**
   - 参数搜索自动化
   - 分析工具增强

## 六、关键决策记录

### Decision 1: 目录结构
- **决策**: 测试放在项目内 `tests/`，与 `ctgomartini/` 平级
- **理由**: 符合 pytest 惯例，方便开发和调试
- **时间**: 2026-03-13

### Decision 2: 测试分层
- **决策**: 分为 unit/integration/e2e 三层
- **理由**: 明确测试范围，unit 快速，e2e 完整但慢
- **时间**: 2026-03-13

### Decision 3: Vermouth 升级策略
- **决策**: 保留 OV+rCSU 为主，martinize2 内部生成为选项
- **理由**: 向后兼容，验证新功能后再切换默认
- **时间**: 2026-03-13

## 七、快速命令

```bash
# 运行测试
pytest tests/                    # 全部
pytest tests/unit/              # 单元测试
pytest tests/integration/       # 集成测试
pytest tests/e2e/               # 端到端测试

# 代码风格（未安装）
# black ctgomartini/
# ruff check ctgomartini/

# 安装
pip install -e ".[dev]"
```

## 八、联系人

- **Author**: Song Yang
- **Email**: yangsong2015@pku.edu.cn
- **Repository**: https://github.com/ComputBiophys/CTGoMartini
