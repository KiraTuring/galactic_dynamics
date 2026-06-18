# 西湖大学计算集群访问配置

## 网络拓扑

```
本地（Windows WSL2）
  └── EasyConnect（Docker 容器）
       ├── noVNC 管理界面: http://localhost:8080（密码: opencode）
       └── SOCKS5 代理: localhost:1080
              └── SSH → 登录节点（10.28.1.66:10002, 用户: wanght245001）
                       ├── sbatch/squeue → 计算节点
                       └── ProxyJump → 计算节点（交互式运行）
```

## 前置条件

- Windows 上安装 Docker Desktop，并在 Settings → Resources → WSL Integration 中启用 Ubuntu
- WSL2 Ubuntu 24.04，Docker Desktop 的 WSL 集成已开启

## 连接步骤

### 1. 启动 VPN

```bash
cd /home/haitong/work/galactic_dynamics/easyconnect
docker compose up -d
```

然后用浏览器打开 http://localhost:8080，输入 VNC 密码 `opencode`，在 EasyConnect 登录界面输入 VPN 账号密码。登录成功后保持容器运行。

### 2. 测试 SSH

```bash
# 通过 SOCKS5 代理连接
ssh -o ProxyCommand='nc -X 5 -x localhost:1080 %h %p' \
    -o StrictHostKeyChecking=accept-new \
    -p 10002 wanght245001@10.28.1.66
```

或使用预配置的 SSH 别名（需将 `easyconnect/ssh_config` 复制到 `~/.ssh/config`）：

```bash
ssh galaxy-login
```

### 3. SSH 配置说明

定义在 `easyconnect/ssh_config` 中的主机别名：

| 别名 | 用途 | 备注 |
|---|---|---|
| `galaxy-login` | 登录节点 | 通过 SOCKS5 代理连接，ControlMaster 持久 4h |
| `galaxy-compute` | 计算节点 | 通过登录节点 ProxyJump，需填写主机名 |

## 集群信息

| 项目 | 内容 |
|---|---|
| 登录节点 | 10.28.1.66:10002 |
| 用户 | wanght245001 |
| 家目录 | `/share/home/maoshudeLab/wanght245001` |
| 项目目录 | `/share/home/maoshudeLab/wanght245001/galactic_dynamics` |
| Python | miniconda3, Python 3.12.9 |
| Slurm 路径 | `/soft/slurm/bin/` |
| 认证方式 | 仅密码（无密钥认证） |

### 分区信息

| 分区 | 用途 | 节点数（空闲） | 备注 |
|---|---|---|---|
| `test` | CPU 通用 | 7 idle | AMD EPYC，默认分区 |
| `test-intel` | CPU 通用 | 16 idle | Intel |
| `l40s` | GPU | 1 idle | NVIDIA L40S |
| `a800` | GPU | 0 idle | NVIDIA A800 |
| `h100` | GPU | 0 idle | NVIDIA H100 |

## Slurm 用法

```bash
# 注意需要用全路径
export PATH=/soft/slurm/bin:$PATH

# 提交作业
sbatch job.sh

# 查看队列
squeue -u wanght245001

# 查看分区
sinfo
```

## 密码自动化

由于集群不支持 SSH 密钥认证，`easyconnect/ssh_helper.py` 提供了自动输密码功能：

```bash
# 只读探查（推荐日常用）
python3 easyconnect/ssh_r.py 'hostname; pwd'
python3 easyconnect/ssh_r.py 'ls /some/dir'
python3 easyconnect/ssh_r.py 'squeue -u wanght245001'

# 安全命令（不需 --exec）
python3 easyconnect/ssh_helper.py 'hostname; pwd'

# 写入命令（需要 --exec）
python3 easyconnect/ssh_helper.py --exec 'mkdir /share/home/.../newdir'

# 自定义超时
python3 easyconnect/ssh_helper.py -t 600 'du -sh /large/dir'
```

**权限模型**：`ssh_r.py` 只允许只读命令（ls, du, cat 等），写入命令硬拒绝。`ssh_helper.py` 允许只读命令，写入命令（rm, mv, mkdir 等）需加 `--exec`（`-x`）。

密码硬编码在脚本中（`SSH_PASSWORD` 环境变量可覆盖）。

## 传输文件

```bash
# 上传（rsync，保留符号链接）
python3 easyconnect/ssh_helper.py --rsync ./local_dir /share/home/.../target/

# scp 方式
scp -o ProxyCommand='nc -X 5 -x localhost:1080 %h %p' \
    -P 10002 local_file wanght245001@10.28.1.66:~/target/
```

## 停止 VPN

```bash
cd /home/haitong/work/galactic_dynamics/easyconnect
docker compose down
```

## 网络故障排查

- Docker Hub 被墙：已配置 DaoCloud 镜像源（`registry-mirrors: ["https://docker.m.daocloud.io"]`）
- EasyConnect 版本不对：尝试 `hagb/docker-easyconnect:7.6.3` 或 `7.6.7`（当前用 7.6.7）
- 需要验证码：CLI 版不支持随机码认证，必须用 VNC/noVNC 图形版
