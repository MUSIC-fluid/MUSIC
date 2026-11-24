# (2+1)D ideal KT prototype for the GPU port

The ideal, boost-invariant KT update of MUSIC (both `face_reconstruction` modes, Heun RK2, the velocity root solve, conformal EoS with n_B). It is written once in `kt2d.h` and compiled as C++/OpenMP, CUDA and Metal. In double it reproduces MUSIC to 5e-9 on the Gubser integration test. See [doc/gpu_porting_investigation.md](../../doc/gpu_porting_investigation.md) for the results and the porting plan.

| file | |
|---|---|
| `kt2d.h` | the kernels (`faceFlux`, `cellUpdate`, `reconstruct`, `maxSpeed`), portable across the three compilers |
| `kt2d_main.cpp` | driver and CPU backend; Gubser initial state, comparison with MUSIC output |
| `kt2d_metal.mm` | Metal backend, MSL compiled at run time from `kt2d.h` |
| `kt2d_cuda.cu` | CUDA backend (`make cuda` with nvcc; `make cuda-syntax` type-checks it with clang) |

```
make && ./kt2d --backend metal --mode 1 --dtau 0.005 --store-q 1
```

The options are listed at the top of `kt2d_main.cpp`.
