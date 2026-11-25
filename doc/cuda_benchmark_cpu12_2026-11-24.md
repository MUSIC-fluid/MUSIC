# CUDA benchmark against the 12-thread CPU backend

On November 24, 2026, the updated MUSIC CUDA build ran **6.96–10.48× faster in FP32** and **1.81–2.03× faster in FP64** than the 12-thread CPU backend on the GTX 1650 laptop. The comparison covers a lumpy TATB event in three hydrodynamic configurations with diagnostics enabled.

## Setup

- **Hardware.** Intel Core i7-9750H with six physical cores and 12 hardware threads, and an NVIDIA GeForce GTX 1650 with 4 GB of device memory.
- **Build.** Release build targeting `sm_75`, using CUDA 13.5 and GCC 15.2. The working tree was based on `e427190` with the CUDA face-kernel EoS specialization and the FP32 hotQCD viscous-cell specialization. All backends used the same executable and shared library.
- **Workload.** The same lumpy TATB initial conditions with baryon density, on a 121×121×37 grid spanning 16 fm × 16 fm × 10 units of spacetime rapidity. Each run completed 61 steps, numbered 0 through 60, at `Delta_Tau = 0.01 fm/c`. The viscous cases used shear with η/s = 0.12, bulk viscosity, and second-order terms. The runs used hydro-only mode with conservation and other standard diagnostics enabled. Freeze-out and full evolution output were disabled.
- **Runtime settings.** The CPU backend used `gpu_backend 0`. CUDA used `gpu_backend 3` with `gpu_float 0` for FP64 and `gpu_float 1` for FP32. The logs confirmed 12 OpenMP threads in every run. CUDA used a requested block size of 128 with kernel profiling disabled.
- **Timing.** Wall time covers the entire MUSIC process, including initialization and diagnostics. There were two timed runs per configuration and 18 timed runs in total. One CUDA warm-up preceded the measurements. Runs were sequential, with workload and backend order reversed in the second pass. Speedup is the CPU median divided by the corresponding CUDA median.

## Wall times

| Workload | Backend | Run 1 (s) | Run 2 (s) | Median (s) | Speedup vs CPU |
|---|---|---:|---:|---:|---:|
| Ideal, EoS 9 | CPU, 12 threads | 22.2719 | 24.7865 | 23.529 | 1.00× |
| Ideal, EoS 9 | CUDA FP64 | 10.5116 | 12.6199 | 11.566 | 2.03× |
| Ideal, EoS 9 | CUDA FP32 | 2.2718 | 2.2206 | 2.246 | 10.48× |
| Viscous, EoS 9 | CPU, 12 threads | 39.4785 | 40.8571 | 40.168 | 1.00× |
| Viscous, EoS 9 | CUDA FP64 | 23.0943 | 21.2399 | 22.167 | 1.81× |
| Viscous, EoS 9 | CUDA FP32 | 4.4328 | 3.9838 | 4.208 | 9.54× |
| Viscous, EoS 20 | CPU, 12 threads | 66.4113 | 66.9917 | 66.701 | 1.00× |
| Viscous, EoS 20 | CUDA FP64 | 31.8966 | 34.7153 | 33.306 | 2.00× |
| Viscous, EoS 20 | CUDA FP32 | 10.0671 | 9.1047 | 9.586 | 6.96× |

## Numerical agreement

Each configuration reproduced its diagnostic files byte for byte across both repeats. The table shows the maximum difference between CUDA and the CPU reference over all 61 conservation records and both CUDA repeats. Energy and net baryon number are normalized by their respective initial CPU values. Zeros indicate agreement to the precision written in the output.

| Workload | CUDA precision | Max `abs(ΔE)/E₀` | Max `abs(ΔN_B)/N_B₀` |
|---|---|---:|---:|
| Ideal, EoS 9 | FP64 | 0 | 0 |
| Ideal, EoS 9 | FP32 | 6.44e-6 | 6.27e-8 |
| Viscous, EoS 9 | FP64 | 0 | 0 |
| Viscous, EoS 9 | FP32 | 6.49e-6 | 6.20e-8 |
| Viscous, EoS 20 | FP64 | 1.41e-10 | 0 |
| Viscous, EoS 20 | FP32 | 6.28e-6 | 5.96e-8 |

The local build directory contains the [raw run records](../build_cuda/benchmark/cpu12_update/runs.jsonl), [executable and input checksums](../build_cuda/benchmark/cpu12_update/build.txt), and [benchmark driver](../build_cuda/benchmark/cpu12_update/compare.py).
