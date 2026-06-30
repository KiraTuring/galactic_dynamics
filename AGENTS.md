# Galactic Dynamics Project

## Overview
Dynamical modeling of galaxies (JAM, AxiSchw, Dynamite).

## Repository Structure

### 本地（本机 WSL2）
```
galactic_dynamics/
├── Axi_Schwarzschild/      # 轴对称 Schwarzschild（含编译好的 Fortran）
├── Trischwarzpy/           # Dynamite wrapper（Dynamite 核心仅在集群）
│   ├── scripts/            # run_bo_dynamite.py, analyze_results.py
│   ├── mod_dyn/            # BayesOpt, 配置, 后处理
│   └── trischwpy/          # Dynamite 包装层 + 三轴数学工具
├── JAM/                    # Jeans Anisotropic MGE 建模
├── data/
│   ├── raw/                # 原始观测数据（git-ignored）
│   └── processed/          # 预处理输出（git-tracked）
├── easyconnect/            # VPN + SSH 工具
├── results/                # 本地分析结果
└── .agents/skills/         # OpenCode skills

### 集群（西湖大学 HPC）
通过 EasyConnect VPN 访问。本地仓库的代码副本通过 `git push/pull` 同步。
```
galactic_dynamics/
├── templates/              # 各星系数据模板（MGE, kinematics）
├── dyn_config/             # Dynamite 模型配置 YAML
├── dyn_models/             # Dynamite 模型输出
├── galaxy_models/          # AxiSchw 模型输出
├── galaxy_data/            # 原始 FITS 数据
├── dynamite/               # Dynamite 源代码（安装版在 schw env 里）
├── _archive/               # 归档数据
├── config/                 # AxiSchw 旧配置
├── JAM/                    # JAM 项目副本
└── notebooks/              # 分析笔记本
```

## Conventions
- Galaxy names: uppercase NGC format (e.g. NGC4621)
- Mock galaxies: append "-mock", "-mock2", etc.
- Custom kinematics: append "-{tag}" to distinguish non-OASIS source (e.g. `-h6` for pPXF h6 fits)
  - no suffix = original OASIS kinematics
  - `-h6` = custom pPXF extraction (moments=6)
- Model configs: `configs/{galaxy}/{model}.yaml`
- Results: `results/JAM/{galaxy}/{model}/`
- Processed data: `data/processed/{galaxy}/` or `data/processed/{galaxy}-{tag}/`

## Data Rules
- Never read FITS files into context — log paths only
- `data/raw/` and `results/` are git-ignored; `data/processed/` IS tracked
- **Never overwrite existing results without explicit user approval.**
  Before writing to `results/` or `data/processed/`, check if data already exists.
  If so, ask the user whether to overwrite, skip, or use a new directory name.
  Scientific data must be traceable and reproducible — accidental deletion
  of hours-long Schwarzschild runs is unacceptable.
- Dynamite 数据流: 本地 `data/processed/` → 上传 → 集群 `templates/` → `dyn_config/xxx.yaml` → `dyn_models/<name>/`

## Axi_Schwarzschild Pipeline
Axisymmetric Schwarzschild orbit-superposition dynamical modeling.
- **Entry point**: `run.py`
- **Core library**: `schwarz.py`
- **Grid generation**: `grid_generate.py` + `gridconfig.yaml`
- **Iterative run**: `iter_generate.py` / `iter_proc.py`
- **Parallel**: `pool_generate.py`
- Requires a compiled `Schwarzschild/` Fortran library.

## JAM — jampy Versions

JAM-fit maintains two branches for different jampy versions:

| Branch | jampy | Notes |
|--------|-------|-------|
| `JAM:` `dev` | 8.1.4 | Legacy MGE API, `logistic=True` keyword |
| `JAM:` `dev-jampy9` | 9.0.0+ | Tuple MGE API, callable `beta`, no `align` |

Switch environment + branch together:

```bash
conda activate schw  && cd JAM && git checkout dev        # jampy 8.1.4
conda activate schw9 && cd JAM && git checkout dev-jampy9  # jampy 9.0.0
```

Full migration notes: `JAM/jampy9-notes.md`.

## Trischwarzpy (Dynamite) Pipeline
Triaxial Schwarzschild orbit-superposition dynamical modeling. Dynamite 核心库仅安装在集群 `schw` conda env 中，本机没有。
- **入口**: `python scripts/run_bo_dynamite.py <config> [-r]`
- **配置**: Dynamite 原生 YAML (`system_components` 格式)
- **并行**: 单节点 multiprocessing（`ncpus` pool workers + `ncpus_weights` NNLS semaphore）
- **恢复**: `-r`/`--resume` 从已有输出继续
- **Git**: 独立仓库，分支 `modified_dyn`，remote `dsimon45/Trischwarzpy.git`
- **关键陷阱**: NNLS 内存爆炸（`ncpus_weights=4` 而非 8）；YAML LaTeX 花括号冲突；Pool 死锁
- 完整文档见 `.agents/skills/dyn-workflow/SKILL.md` 和 `Trischwarzpy/AGENTS.md`

## Environment
- Python 3.11+
- Key dependencies: jampy, dynesty, cmaes, adamet, mgefit, plotbin, astropy
- 集群 conda env: `schw`（含 dynamite, pathos, skopt, mpi4py 等）
- 注意: 集群 `python3` 是 base env，Dynamite 必须用 `schw` env
  ```
  /share/.../miniconda3/envs/schw/bin/python3
  ```

## 计算集群（西湖大学）

通过 EasyConnect VPN + SOCKS5 代理访问。完整文档见 `easyconnect/CLUSTER.md`。

### 快速参考

```bash
# 启动 VPN（浏览器 http://localhost:8080 登录）
cd easyconnect && docker compose up -d

# SSH 登录（密码: R3w44CWWc*GATbk，或用 ssh_helper.py 自动输入）
ssh -o ProxyCommand='nc -X 5 -x localhost:1080 %h %p' -p 10002 wanght245001@10.28.1.66

# 执行命令（自动输密码）
python3 easyconnect/ssh_helper.py '/soft/slurm/bin/squeue -u wanght245001'

# 提交作业
python3 easyconnect/ssh_helper.py 'export PATH=/soft/slurm/bin:$PATH && sbatch job.sh'

# 写命令需加 --exec/-x
python3 easyconnect/ssh_helper.py --exec 'sbatch job.sh'

# 上传/下载文件
python3 easyconnect/ssh_helper.py --push local_file cluster_path
python3 easyconnect/ssh_helper.py --pull cluster_path local_path
python3 easyconnect/ssh_helper.py --download cluster_file local_path
```

### 主机别名（`~/.ssh/config`）
- `galaxy-login` — 登录节点（SOCKS5 代理）
- `galaxy-compute` — 计算节点（ProxyJump 通过登录节点）

### 集群路径
- 项目目录: `/share/home/maoshudeLab/wanght245001/galactic_dynamics/`
- Slurm: `/soft/slurm/bin/`
- Python: `/share/home/maoshudeLab/wanght245001/miniconda3/bin/python3`
- Python (schw env): `/share/home/maoshudeLab/wanght245001/miniconda3/envs/schw/bin/python3`

### 关键文件
| 文件 | 用途 |
|---|---|
| `easyconnect/docker-compose.yml` | EasyConnect 容器配置 |
| `easyconnect/ssh_helper.py` | 自动输密码 SSH |
| `easyconnect/ssh_config` | SSH 配置模板 |
| `easyconnect/CLUSTER.md` | 完整集群文档 |
| `opencode.jsonc` | 项目配置（含 cluster 元信息） |

## Skills
- `galaxy-data-prep` — 数据预处理（OASIS/SAURON kinematics, MGE, PSF, custom pPXF）
- `galaxy-jam` — JAM config 创建、运行、结果管理
- `galaxy-troubleshooting` — 管线常见问题排查
- `mge-fitting` — 从 HST 图像拟合 MGE 光度模型（Cappellari 2002, mgefit）
- `schw-workflow` — AxiSchw 完整工作流（配置 → 提交 → 监控 → 分析）
- `dyn-workflow` — Dynamite 三轴完整工作流（配置 → 提交 → 恢复 → 后处理）
