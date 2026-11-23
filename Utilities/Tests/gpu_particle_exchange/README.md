# GPU 粒子跨网格交换回归

本组检查 G 路径中液滴的归属判断、删除、压缩、打包和接收；不是完整的粒子模型验证套件。
输入使用 VLES、无随机扰动和开放边界，关闭液滴传热传质，以便单独检验迁移和守恒。

## 实现范围

入口为 `G_EXCH::GPU_MESH_EXCHANGE_PARTICLES`，实现位于 `Source/g_part.f90`。

- GPU 进行目标网格分类、寿命/半径/离域删除、周期回绕、打包、压缩和接收。
- 同 rank 的网格直接读取设备发送缓冲；跨 rank 使用既有 CUDA-aware MPI 包装传递设备缓冲。
  整型 TAG、类别和元数据使用整型缓冲，不转换为浮点数。
- CPU 只调度 kernel/MPI、处理通信计数和缓冲容量；不逐粒子计算迁移。
  缓冲按容量复用，必要时扩容并在 GPU 内复制。
- 独立 MPI communicator 隔离粒子消息；先发布全部接收，再发送，支持一个 rank 管理多个网格。
- 设备粒子数组是计算状态；CPU 粒子数组只是输出镜像及新注入批次的暂存区。
  不再通过每步整批 D2H、CPU `REMOVE_PARTICLES`、CPU 交换、GPU 重新初始化完成迁移。
- `NUMBER OF PARTICLES` 使用交换后已有的计数，不依赖输出镜像，也不要求整批粒子回传。
- C 和 GC 保留 CPU 粒子计算/交换；GC 测试隐藏 CUDA 设备，检查它不依赖 GPU 初始化。

## 算例

| 输入 | 检查内容 | 最终各网格粒子数 |
| --- | --- | --- |
| `one_way.fds` | 单向迁入初始为空的网格，接收容量扩展 | 0, 32 |
| `bidirectional.fds` | 三网格、双向、两种液滴类别 | 21, 43, 22 |
| `reinsertion.fds` | 迁移后继续注入，粒子输出不逐步进行 | 4, 28 |
| `removal.fds` | 寿命与离域删除，允许最终空网格 | 0, 0 |
| `periodic.fds` | 经过另一网格后从周期对侧返回 | 32, 0 |
| `linked.fds` | INIT 关联的 DEVC 随粒子跨网格继续跟踪 | 0, 1 |
| `many_particles.fds` | 649 个粒子，多 CUDA block；512 个迁入细网格，137 个留下 | 137, 512 |
| `multidest.fds` | 四网格；一个源同时向三个目标发送，含同一 rank 的多个目标 | 0, 40, 41, 42 |
| `count_devices.fds` | 每步读取网格粒子数，粒子镜像仅间隔输出；与 C 对照 | 4, 28 |

`check_particles.py` 读取真实 `.prt5` 二进制记录，按 TAG 而非数组顺序比较。
检查粒子无重复、无非预期丢失，类别不变，位置落在所属网格，质量/权重守恒，
迁移和删除达到预期结果；随后与 C 的位置和物理量比较。
`count_devices` 还逐行比较 DEVC 计数，检查输出镜像未刷新时的计数。
程序退出码不是唯一通过条件：日志必须含正常结束标记且不含 `ERROR`。

`.prt5` 中数值为单精度；这里的“输出差为零”不等于内部双精度位级一致。
这些小算例不用于评价多 GPU 加速比，也不能替代有蒸发、辐射或湍流时的物理验证。

## 运行

需要已编译的 GFDS、支持 CUDA-aware MPI 的运行环境以及 Python 3 标准库，无额外 Python 包。
务必等构建/链接完全结束后再执行。脚本严格串行运行算例，不同时占用 GPU。

```bash
cmake --build run -j4
bash Utilities/Tests/gpu_particle_exchange/run.sh "$PWD/run/fds"
```

默认使用新临时目录；也可指定**尚不存在**的结果目录及已有 C 基准根目录：

```bash
bash Utilities/Tests/gpu_particle_exchange/run.sh \
  "$PWD/run/fds" /path/to/new-results /path/to/previous-results
```

每例执行 C、GC、G1、G：G1 是一个 rank 管理多个网格；G 通常为两个 rank，
`bidirectional` 为三个 rank。GC 设置 `CUDA_VISIBLE_DEVICES=-1`。
只有输入逐字一致且检查器通过时才复用 C 基准，并以符号链接保留原结果；
新增或变更输入会单独生成 C 基准。已存在的目标结果目录会被拒绝，不覆盖旧输出。

单独检查已有结果：

```bash
python3 Utilities/Tests/gpu_particle_exchange/check_particles.py \
  /path/to/results/many_particles_G --case many_particles \
  --reference /path/to/results/many_particles_C
```

## 本机验证证据（2026-11-23）

结果根目录为 `run/test/g_part/device_exchange_20261123/`（运行产物不纳入 Git）。

- `final_checked/`：补入当前计数接口后，九例 C/GC/G1/G 全部检查通过；前八例复用
  输入相同且已检查的 C 基准，仅新增 `count_devices` 的 C 运行。
  `count_devices` 每 0.002 s 的 DEVC 计数与 C 一致，粒子文件间隔为 0.05 s，
  检查覆盖两次粒子输出之间的注入和迁移。
- `release_checked/`：前八例的 C/GC/G1/G 检查通过。大部分 C 基准复用已有结果。
  与 C 比较，输出位置、质量、权重和温度的最大差为零；重复注入例的最大 U 差约
  `4.77e-7 m/s`，多目标例约 `3.26e-14 m/s`。
- `memcheck_release/`、`initcheck_release/`：`many_particles` 的设备内存检查均为 0 errors，
  粒子记录检查同时通过。检查采用下述单 rank MPI 配置；不代表完成跨 rank sanitizer 检查。
- `residency_release/particle_exchange.nsys-rep` 和 `.sqlite`：在 512 个迁移、137 个留下的
  一次同 rank 交换窗口内，从分类到接收完成仅记录 16 字节计数 D2H、8 字节偏移 H2D，
  没有整批粒子物理量回传。这是该窗口的证据，不代表整个求解过程完全没有 D2H，
  也不能证明 MPI 库底层没有内部 staging。
- `LN02_4_release/`：5 s、2500 步短回归正常结束，后段每百步约 1.06–1.10 s。
  对同输入的已有 C 结果，14 条 LINE 曲线相对 L2 差约 1.15%–5.80%，
  最低相关系数约 0.9973。不是完整 50 s 实验验证，差异不能直接归为浮点舍入。
  正常结束时仍有 divide-by-zero/underflow/inexact 浮点标志，本项尚未追查其来源。

内存检查应在独立的新算例目录进行，避免覆盖结果，例如在复制了 `many_particles.fds` 的目录运行：

```bash
mpirun --mca osc pt2pt --mca pml ob1 --mca btl self,vader -np 1 \
  compute-sanitizer --tool memcheck --error-exitcode 86 --log-file memcheck.log \
  /absolute/path/to/fds many_particles.fds --g
```

将 `memcheck` 换为 `initcheck` 并更换日志名即可做初始化检查。
本机默认 UCX 与 sanitizer 组合曾在 `MPI_Init` 中报无 CUDA context，故使用上述配置；
该失败不能作为通过结果。普通多 rank 运行仍有 `cuMemHostRegister ... failed` 环境警告。

## 边界与未完成工作

- 实测为一张 RTX 4060 Laptop GPU 上的多个 MPI rank，未验证真实多 GPU、多节点或 GPUDirect 性能。
- 连续喷淋/体积粒子生成 `INSERT_ALL_PARTICLES` 仍在 CPU，只上传新批次。
- PDPA 和部分粒子 DEVC 仍使用 CPU 镜像，可能每步回传粒子；不能把它们统称为低频文件输出。
- CTRL 的 CPU `UPDATE_CONTROLS`、HVAC 的 CPU `HVAC_CALC` 不属于本次交换改造。
- 保留已有液滴模型限制：固壁相互作用、破碎、湍流弥散、路径/风管粒子、柱坐标、
  `PARTICLE_CFL` 等未实现的 G 分支仍拒绝运行，不静默退回 CPU 算法。
- 达到 `MAXIMUM_PARTICLES` 后“删除最旧粒子”的 G 策略尚未实现，当前明确报错；C/GC 不受此限制影响。
- restart 写出前已同步粒子镜像，但完整的 restart 读写循环尚无专项回归。

本组通过只能说明当前支持范围内的 GPU 粒子迁移经过上述检查，不代表所有粒子/CTRL/HVAC 算法已全 GPU。
