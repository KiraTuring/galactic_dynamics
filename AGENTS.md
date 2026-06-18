# Galactic Dynamics Project

## Overview
Dynamical modeling of galaxies (JAM, AxiSchw, future: Dynamite).

## Repository Structure
```
galactic_dynamics/
├── data/
│   ├── raw/                # Original FITS & raw tables (per galaxy)
│   └── processed/          # MGE, aperture, bins, kinematics (per galaxy)
├── shared_utils/           # Shared preprocessing utilities (future)
├── JAM/
│   ├── jam_fit/            # Core JAM library
│   ├── configs/{galaxy}/   # Model config YAML files
│   ├── scripts/            # Main entry points & utilities
│   ├── notebooks/          # Exploration & visualization
│   └── AGENTS.md           # JAM-specific conventions
├── Axi_Schwarzschild/      # Axisymmetric Schwarzschild orbit-superposition code
├── results/
│   └── JAM/{galaxy}/{model}/
├── .gitignore
└── README.md
```

## Conventions
- Galaxy names: uppercase NGC format (e.g. NGC4621)
- Mock galaxies: append "-mock", "-mock2", etc.
- Model configs: `configs/{galaxy}/{model}.yaml`
- Results: `results/JAM/{galaxy}/{model}/`
- Processed data: `data/processed/{galaxy}/`

## Data Rules
- Never read FITS files into context — log paths only
- `data/raw/` and `results/` are git-ignored
- Dynamite data lives in `data/processed/{galaxy}-dyn/` and `{galaxy}-gh/`

## Axi_Schwarzschild Pipeline
Axisymmetric Schwarzschild orbit-superposition dynamical modeling.
- **Entry point**: `run.py`
- **Core library**: `schwarz.py`
- **Grid generation**: `grid_generate.py` + `gridconfig.yaml`
- **Iterative run**: `iter_generate.py` / `iter_proc.py`
- **Parallel**: `pool_generate.py`
- Requires a compiled `Schwarzschild/` Fortran library.

## Environment
- Python 3.11+
- Key dependencies: jampy, dynesty, cmaes, adamet, mgefit, plotbin, astropy

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
```

### 主机别名（`~/.ssh/config`）
- `galaxy-login` — 登录节点（SOCKS5 代理）
- `galaxy-compute` — 计算节点（ProxyJump 通过登录节点）

### 集群路径
- 项目目录: `/share/home/maoshudeLab/wanght245001/galactic_dynamics/`
- Slurm: `/soft/slurm/bin/`
- Python: `/share/home/maoshudeLab/wanght245001/miniconda3/bin/python3`

### 关键文件
| 文件 | 用途 |
|---|---|
| `easyconnect/docker-compose.yml` | EasyConnect 容器配置 |
| `easyconnect/ssh_helper.py` | 自动输密码 SSH |
| `easyconnect/ssh_config` | SSH 配置模板 |
| `easyconnect/CLUSTER.md` | 完整集群文档 |
| `opencode.jsonc` | 项目配置（含 cluster 元信息） |
