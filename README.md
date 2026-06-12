# FMCW Radar Signal Processing — Multi-Target Simulation

A correct, well-documented MATLAB simulation of FMCW (Frequency-Modulated Continuous Wave) radar range-Doppler processing for multi-target scenarios. Designed as the algorithmic foundation for a future CUDA-accelerated radar signal processing library.

## FMCW Principle

### Waveform

FMCW radar transmits a continuous wave whose frequency linearly increases with time. Each frequency sweep is called a *chirp*:

$$s_{tx}(t) = \exp\{j2\pi[f_c t + \frac{1}{2}S\tau^2]\}$$

where $f_c$ is the carrier frequency, $S = B/T_{chirp}$ is the slope (sweep rate), $\tau = t \bmod T_{chirp}$ is the fast time within a chirp, and $T_{chirp}$ is the chirp duration.

### Fast-Time / Slow-Time Framework

Signal processing is organized on a 2D grid:

| Dimension | Time scale | Samples | What it captures |
|-----------|-----------|---------|-----------------|
| **Fast time** $\tau$ ($N_r$) | $0 \sim T_{chirp}$ | Range bins | Distance (via IF frequency) |
| **Slow time** $n$ ($N_d$) | $0 \sim N_d \cdot T_{chirp}$ | Chirp index | Velocity (via Doppler phase) |

### Range Measurement (Dechirp-on-Receive)

The received signal from a target at delay $\tau_k = 2R_k/c$ is:

$$s_{rx}(t) = \exp\{j2\pi[f_c(t-\tau_k) + \frac{1}{2}S(\tau-\tau_k)^2]\}, \quad \tau \ge \tau_k$$

Mixing (multiplying TX by conjugate of RX) removes the carrier and chirp, leaving the beat signal:

$$s_{if} = s_{tx} \cdot s_{rx}^* = \exp\{j2\pi[\underbrace{f_c\tau_k}_{\text{carrier phase}} + \underbrace{S\tau\tau_k}_{\text{beat frequency}} - \underbrace{\frac{1}{2}S\tau_k^2}_{\text{RVP}}]\}$$

The beat frequency $f_{if} = S \cdot \tau_k = \frac{B}{T_{chirp}} \cdot \frac{2R}{c}$ is proportional to range. A 1D FFT along fast time extracts the range.

### Velocity Measurement (Doppler)

The carrier phase term $2\pi f_c \tau_k$ is constant within a single chirp but varies slowly across chirps as the target moves. A 2D FFT along slow time extracts this phase variation as Doppler frequency:

$$f_d = \frac{2v}{\lambda}, \quad v = \frac{\lambda}{2} f_d$$

### Key Performance Limits

| Metric | Formula | This simulation |
|--------|---------|-----------------|
| Range resolution | $\Delta R = c/(2B)$ | 1.0 m |
| Max unambiguous range (Nyquist) | $R_{max} = \frac{N_r}{2} \cdot \Delta R$ | 512 m |
| Max unambiguous range (time delay) | $R_{max}^{abs} = cT_{chirp}/2$ | 7500 m |
| Velocity resolution | $\Delta v = \lambda/(2 N_d T_{chirp})$ | 0.15 m/s |
| Max unambiguous velocity | $v_{max} = \lambda/(4 T_{chirp})$ | 19.5 m/s |

> **Note on max range:** The Nyquist limit ($R_{max} = 512$ m) is the stricter constraint. Beyond this distance, the IF frequency exceeds half the sampling rate, causing aliasing. The time-delay limit ($R_{max}^{abs} = 7500$ m) is the absolute physical maximum — beyond it, the echo arrives after the chirp ends and no samples are valid.

## Features

- ✅ **Physically correct** — fast/slow time separation, stop-and-hop, RX delayed naturally
- ✅ **Multi-target** — configurable positions, velocities, amplitudes
- ✅ **Vectorized** — no nested per-sample loops
- ✅ **Controllable noise** — beat-level SNR with 2D FFT processing gain
- ✅ **Windowing** — Hanning window to suppress FFT sidelobes (-32 dB vs -13 dB)
- ✅ **CA-CFAR detection** — 2D sliding window in power domain
- ✅ **DBSCAN clustering** — normalized isotropic coordinates, power-weighted centroids
- ✅ **IEEE-style visualization** — 3D RDM, plan view, range/Doppler profiles, CFAR overlay

## File Structure

```
FMCW/
├── fmcw_simulation.m       # Main simulation script
├── README.md               # This file
└── LICENSE                  # MIT License
```

## Usage

Run directly in MATLAB:

```matlab
>> fmcw_simulation
```

The script will:

1. Parse radar parameters and display them
2. Define 5 targets with varying range/velocity
3. Generate the beat signal (2D matrix: $N_r \times N_d$)
4. Perform 2D FFT with Hanning windowing → Range-Doppler Map
5. Run CA-CFAR detection
6. Apply DBSCAN clustering to aggregate detections
7. Display 7 figures: 3D RDM, 2D plan view, range profile, Doppler profile, and CFAR results

### Configuration

Edit the parameters at the top of `fmcw_simulation.m`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `fc` | 77e9 | Carrier frequency (Hz) |
| `B` | 150e6 | Sweep bandwidth (Hz) |
| `T_chirp` | 50e-6 | Chirp duration (s) |
| `Nd` | 256 | Number of chirps per frame |
| `Nr` | 1024 | Samples per chirp |
| `SNR_beat_dB` | -10 | Beat-level SNR (dB) |
| `offset` | 13 | CFAR threshold offset (dB) |

## Radar Parameters (Default)

| Parameter | Value |
|-----------|-------|
| Carrier frequency | 77 GHz |
| Bandwidth | 150 MHz |
| Chirp duration | 50 μs |
| Range resolution | 1.0 m |
| Max range (Nyquist) | 512 m |
| Velocity resolution | 0.15 m/s |
| Max velocity | 19.5 m/s |

## Targets (Default)

| Target | Range | Velocity | Amp |
|--------|-------|----------|-----|
| T1 | 30 m | 0 m/s | 1.0 |
| T2 | 50 m | +5 m/s | 0.8 |
| T3 | 75 m | -10 m/s | 1.2 |
| T4 | 120 m | 0 m/s | 0.5 |
| T5 | 155 m | +8 m/s | 0.6 |

## License

MIT
