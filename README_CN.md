# FMCW 雷达信号处理 — 多目标仿真

一个正确的、注释完善的 FMCW（调频连续波）雷达多目标 Range-Doppler 处理 MATLAB 仿真。本仿真是后续 CUDA 加速雷达信号处理算法库的信号模型基础。

## FMCW 原理

### 波形

FMCW 雷达发射频率随时间线性上升的连续波信号。每一次频率扫描称为一个 chirp（线性调频脉冲）：

$$s_{tx}(t) = \exp\{j2\pi[f_c t + \frac{1}{2}S\tau^2]\}$$

其中 $f_c$ 为载频，$S = B/T_{chirp}$ 为调频斜率，$\tau = t \bmod T_{chirp}$ 为 chirp 内的快时间，$T_{chirp}$ 为 chirp 时长。

### 快时间 / 慢时间框架

信号处理在二维网格上组织：

| 维度 | 时间尺度 | 采样数 | 捕获的信息 |
|------|---------|--------|-----------|
| **快时间** $\tau$ | $0 \sim T_{chirp}$ | 距离 bin $N_r$ | 距离（通过中频频率） |
| **慢时间** $n$ | $0 \sim N_d T_{chirp}$ | chirp 数 $N_d$ | 速度（通过多普勒相位） |

### 测距原理（Dechirp 接收）

时延为 $\tau_k = 2R_k/c$ 的目标，其接收信号为：

$$s_{rx}(t) = \exp\{j2\pi[f_c(t-\tau_k) + \frac{1}{2}S(\tau-\tau_k)^2]\}, \quad \tau \ge \tau_k$$

发射与接收共轭相乘（复混频），去掉载波和 chirp 二次项，得到差频信号：

$$s_{if} = s_{tx} \cdot s_{rx}^* = \exp\{j2\pi[\underbrace{f_c\tau_k}_{\text{载波相位}} + \underbrace{S\tau\tau_k}_{\text{差频项}} - \underbrace{\frac{1}{2}S\tau_k^2}_{\text{RVP}}]\}$$

差频频率 $f_{if} = S \cdot \tau_k = \dfrac{B}{T_{chirp}} \cdot \dfrac{2R}{c}$ 与距离成正比。沿快时间做 FFT 即可提取距离。

### 测速原理（多普勒）

载波相位项 $2\pi f_c \tau_k$ 在单 chirp 内为常数，但跨 chirp 时会因目标运动而线性变化。沿慢时间做 FFT 提取该相位变化，得到多普勒频率：

$$f_d = \frac{2v}{\lambda}, \quad v = \frac{\lambda}{2} f_d$$

### 关键性能极限

| 指标 | 公式 | 本仿真参数 |
|------|------|-----------|
| 距离分辨率 | $\Delta R = c/(2B)$ | 1.0 m |
| 最大无模糊距离（Nyquist） | $R_{max} = \dfrac{N_r}{2} \cdot \Delta R$ | 512 m |
| 最大无模糊距离（时延） | $R_{max}^{abs} = cT_{chirp}/2$ | 7500 m |
| 速度分辨率 | $\Delta v = \lambda/(2 N_d T_{chirp})$ | 0.15 m/s |
| 最大无模糊速度 | $v_{max} = \lambda/(4 T_{chirp})$ | 19.5 m/s |

> **关于最大探测距离的说明：** Nyquist 极限（512 m）是更严格的约束。超过此距离时，中频频率超过采样率的一半（$F_s/2$），产生频谱混叠。时延极限（7500 m）是物理极限——超过此距离时，回波在 chirp 结束后才到达，有效采样点数为 0。

## 功能特性

- ✅ **物理正确** — 快-慢时间分离、stop-and-hop 近似、RX 信号自动适配时延
- ✅ **多目标** — 可配置位置、速度、幅度
- ✅ **向量化加速** — 无逐采样点嵌套循环
- ✅ **可控噪声** — beat 层面定义 SNR，含 2D FFT 处理增益
- ✅ **加窗** — Hanning 窗将旁瓣从 -13 dB 压制到 -32 dB
- ✅ **CA-CFAR 检测** — 二维滑动窗，功率域处理
- ✅ **DBSCAN 聚类** — 归一化各向同性坐标 + 功率加权质心
- ✅ **IEEE 风格可视化** — 3D RDM、平面图、距离/多普勒剖面、CFAR 叠加

## 文件结构

```
FMCW/
├── fmcw_simulation.m       # 主仿真脚本
├── README.md               # 英文说明
├── README_CN.md            # 中文说明
├── LICENSE                 # MIT License
└── .gitignore              # Git 忽略规则
```

## 使用方法

直接在 MATLAB 中运行：

```matlab
>> fmcw_simulation
```

脚本自动完成以下步骤：

1. 显示雷达参数
2. 定义 5 个目标（距离/速度/幅度）
3. 生成差频信号（$N_r \times N_d$ 二维矩阵）
4. Hanning 加窗 + 2D FFT → Range-Doppler Map
5. CA-CFAR 检测
6. DBSCAN 聚类聚合检测点
7. 输出 7 张图：3D RDM、2D 平面图、距离剖面、多普勒剖面、CFAR 结果

### 参数配置

修改 `fmcw_simulation.m` 顶部的参数：

| 参数 | 默认值 | 说明 |
|------|--------|------|
| `fc` | 77e9 | 载频 (Hz) |
| `B` | 150e6 | 带宽 (Hz) |
| `T_chirp` | 50e-6 | Chirp 时长 (s) |
| `Nd` | 256 | 每帧 chirp 数 |
| `Nr` | 1024 | 每 chirp 采样数 |
| `SNR_beat_dB` | -10 | Beat 层面信噪比 (dB) |
| `offset` | 13 | CFAR 门限偏置 (dB) |

## 默认雷达参数

| 参数 | 值 |
|------|----|
| 载频 | 77 GHz |
| 带宽 | 150 MHz |
| Chirp 时长 | 50 μs |
| 距离分辨率 | 1.0 m |
| 最大距离（Nyquist） | 512 m |
| 速度分辨率 | 0.15 m/s |
| 最大速度 | 19.5 m/s |

## 默认目标

| 目标 | 距离 | 速度 | 幅度 |
|------|------|------|------|
| T1 | 30 m | 0 m/s | 1.0 |
| T2 | 50 m | +5 m/s | 0.8 |
| T3 | 75 m | -10 m/s | 1.2 |
| T4 | 120 m | 0 m/s | 0.5 |
| T5 | 155 m | +8 m/s | 0.6 |

## License

MIT
