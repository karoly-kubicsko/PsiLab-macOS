// ΨLab 15.0: optional acceleration of the matrix products (gemm in linalg.h), loaded at run time.
//
//  - GPU offload of large products: NVIDIA through cuBLAS (libcublas + libcudart), AMD through
//    hipBLAS (libhipblas + libamdhip64, ROCm). A product is sent to the GPU when it has at least
//    PSILAB_GPU_MIN_FLOPS floating-point operations (default 2e9, 2 m n k): below that the copies
//    over PCIe cost more than the product. Double precision throughout (the energies need it);
//    data-centre GPUs (NVIDIA A100/H100/B200, AMD Instinct MI200/MI300) have full-rate FP64,
//    consumer GPUs 1/32 or 1/64 of it, where the CPU is often as fast.
//  - An optimised CPU BLAS on Linux and Windows: OpenBLAS, Intel MKL (mkl_rt), BLIS or the
//    system's libblas (cblas_dgemm, also the 64-bit-integer and SciPy-prefixed builds of
//    OpenBLAS). On macOS the backend uses Apple's Accelerate (compiled in), which runs on the
//    matrix units of Apple silicon; Apple GPUs have no FP64, so there is no GPU path on macOS.
//
// Environment variables: PSILAB_GPU = auto (default) | off | cuda | hip; PSILAB_BLAS = auto
// (default) | off | <path of a BLAS library>; PSILAB_GPU_MIN_FLOPS = <number>. Nothing is
// required: without these libraries the portable blocked GEMM runs, multi-threaded.
#pragma once

#include <cstddef>
#include <string>

namespace macqc {

// Tries the accelerated paths for C = alpha op(A) op(B) + beta C (row-major, as gemm()); false
// if none applies (the caller then multiplies itself).
bool accelerated_gemm(bool trans_a, bool trans_b, std::size_t m, std::size_t n, std::size_t k, double alpha,
                      const double* a, std::size_t lda, const double* b, std::size_t ldb, double beta, double* c,
                      std::size_t ldc, int threads);

// What is available: e.g. "NVIDIA GPU (cuBLAS, 1 device) for products >= 2e+09 flop; OpenBLAS
// (libopenblas.so.0)" or "" (nothing loaded).
std::string accelerator_description();

// Counters (for tests and the output file): products done on the GPU and by the CPU BLAS.
struct AcceleratorCounters {
    std::size_t gpu_products = 0, blas_products = 0;
    double gpu_flops = 0.0;
};
AcceleratorCounters accelerator_counters();

}  // namespace macqc
