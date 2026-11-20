// ΨLab 15.0: run-time loaded GPU (cuBLAS, hipBLAS) and CPU BLAS paths of gemm(); see
// macqc/accelerators.h.
#include "macqc/accelerators.h"

#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <string>
#include <vector>

#if !defined(_WIN32)
#include <dlfcn.h>
#endif

namespace macqc {

namespace {

std::string env(const char* name) {
    const char* v = std::getenv(name);
    std::string s = v ? v : "";
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

void* open_library(const std::vector<const char*>& names, std::string* which) {
#if defined(_WIN32)
    (void)names;
    (void)which;
    return nullptr;
#else
    for (const char* n : names) {
        if (void* h = dlopen(n, RTLD_NOW | RTLD_LOCAL)) {
            if (which) *which = n;
            return h;
        }
    }
    return nullptr;
#endif
}

void* symbol(void* lib, const char* name) {
#if defined(_WIN32)
    (void)lib;
    (void)name;
    return nullptr;
#else
    return lib ? dlsym(lib, name) : nullptr;
#endif
}

// ---- CPU BLAS (cblas_dgemm with 32-bit or 64-bit integers)

enum { kRowMajor = 101, kNoTrans = 111, kTrans = 112 };
using cblas32_t = void (*)(int, int, int, int, int, int, double, const double*, int, const double*, int, double,
                           double*, int);
using cblas64_t = void (*)(int, int, int, std::int64_t, std::int64_t, std::int64_t, double, const double*, std::int64_t,
                           const double*, std::int64_t, double, double*, std::int64_t);
using set_threads_t = void (*)(int);

struct CpuBlas {
    cblas32_t g32 = nullptr;
    cblas64_t g64 = nullptr;
    set_threads_t set_threads = nullptr;
    std::string name;
};

CpuBlas load_cpu_blas() {
    CpuBlas b;
#if defined(MACQC_USE_ACCELERATE)
    return b;  // Accelerate is compiled in
#else
    const std::string mode = env("PSILAB_BLAS");
    if (mode == "off" || mode == "no" || mode == "0") return b;
    std::vector<const char*> names;
    const char* explicit_path = std::getenv("PSILAB_BLAS");
    if (explicit_path && *explicit_path && mode != "auto" && mode != "on") names.push_back(explicit_path);
    for (const char* n : {"libopenblas.so.0", "libopenblas.so", "libmkl_rt.so.2", "libmkl_rt.so", "libblis.so.4",
                          "libblis.so", "libopenblas.dylib"})  // not the reference libblas: slower than the portable GEMM
        names.push_back(n);
    std::string which;
    void* lib = open_library(names, &which);
    if (!lib) return b;
    if (auto f = reinterpret_cast<cblas32_t>(symbol(lib, "cblas_dgemm"))) b.g32 = f;
    else if (auto f2 = reinterpret_cast<cblas32_t>(symbol(lib, "scipy_cblas_dgemm"))) b.g32 = f2;
    else if (auto f3 = reinterpret_cast<cblas64_t>(symbol(lib, "cblas_dgemm64_"))) b.g64 = f3;
    else if (auto f4 = reinterpret_cast<cblas64_t>(symbol(lib, "scipy_cblas_dgemm64_"))) b.g64 = f4;
    if (!b.g32 && !b.g64) return CpuBlas{};
    b.set_threads = reinterpret_cast<set_threads_t>(symbol(lib, "openblas_set_num_threads"));
    if (!b.set_threads) b.set_threads = reinterpret_cast<set_threads_t>(symbol(lib, "MKL_Set_Num_Threads"));
    b.name = which;
    return b;
#endif
}

// ---- GPU (cuBLAS or hipBLAS, with the runtime's memory functions)

using create_t = int (*)(void**);
using dgemm_t = int (*)(void*, int, int, int, int, int, const double*, const double*, int, const double*, int,
                        const double*, double*, int);
using malloc_t = int (*)(void**, std::size_t);
using free_t = int (*)(void*);
using memcpy_t = int (*)(void*, const void*, std::size_t, int);
using count_t = int (*)(int*);

struct Gpu {
    bool ok = false;
    std::string kind;   // "cuda" or "hip"
    std::string name;
    int devices = 0;
    void* handle = nullptr;
    dgemm_t dgemm = nullptr;
    malloc_t dmalloc = nullptr;
    free_t dfree = nullptr;
    memcpy_t dmemcpy = nullptr;
    int op_n = 0, op_t = 1, to_device = 1, to_host = 2;
    double min_flops = 2.0e9;
    // device buffers, grown as needed
    void* buf[3] = {nullptr, nullptr, nullptr};
    std::size_t cap[3] = {0, 0, 0};
};

Gpu load_gpu() {
    Gpu g;
    const std::string mode = env("PSILAB_GPU");
    if (mode == "off" || mode == "no" || mode == "0") return g;
    if (const char* f = std::getenv("PSILAB_GPU_MIN_FLOPS")) {
        const double v = std::atof(f);
        if (v > 0.0) g.min_flops = v;
    }
    auto try_cuda = [&]() {
        std::string rt, bl;
        void* cudart = open_library({"libcudart.so", "libcudart.so.13", "libcudart.so.12", "libcudart.so.11.0"}, &rt);
        void* cublas = open_library({"libcublas.so", "libcublas.so.13", "libcublas.so.12", "libcublas.so.11"}, &bl);
        if (!cudart || !cublas) return false;
        auto count = reinterpret_cast<count_t>(symbol(cudart, "cudaGetDeviceCount"));
        auto create = reinterpret_cast<create_t>(symbol(cublas, "cublasCreate_v2"));
        g.dgemm = reinterpret_cast<dgemm_t>(symbol(cublas, "cublasDgemm_v2"));
        g.dmalloc = reinterpret_cast<malloc_t>(symbol(cudart, "cudaMalloc"));
        g.dfree = reinterpret_cast<free_t>(symbol(cudart, "cudaFree"));
        g.dmemcpy = reinterpret_cast<memcpy_t>(symbol(cudart, "cudaMemcpy"));
        if (!count || !create || !g.dgemm || !g.dmalloc || !g.dfree || !g.dmemcpy) return false;
        int n = 0;
        if (count(&n) != 0 || n < 1) return false;
        if (create(&g.handle) != 0) return false;
        g.devices = n;
        g.kind = "cuda";
        g.name = "NVIDIA GPU (cuBLAS, " + bl + ")";
        g.op_n = 0;  // CUBLAS_OP_N
        g.op_t = 1;  // CUBLAS_OP_T
        g.to_device = 1;  // cudaMemcpyHostToDevice
        g.to_host = 2;    // cudaMemcpyDeviceToHost
        return true;
    };
    auto try_hip = [&]() {
        std::string rt, bl;
        void* hiprt = open_library({"libamdhip64.so", "libamdhip64.so.6", "libamdhip64.so.7", "libamdhip64.so.5"}, &rt);
        void* hipblas = open_library({"libhipblas.so", "libhipblas.so.2", "libhipblas.so.3", "libhipblas.so.1"}, &bl);
        if (!hiprt || !hipblas) return false;
        auto count = reinterpret_cast<count_t>(symbol(hiprt, "hipGetDeviceCount"));
        auto create = reinterpret_cast<create_t>(symbol(hipblas, "hipblasCreate"));
        g.dgemm = reinterpret_cast<dgemm_t>(symbol(hipblas, "hipblasDgemm"));
        g.dmalloc = reinterpret_cast<malloc_t>(symbol(hiprt, "hipMalloc"));
        g.dfree = reinterpret_cast<free_t>(symbol(hiprt, "hipFree"));
        g.dmemcpy = reinterpret_cast<memcpy_t>(symbol(hiprt, "hipMemcpy"));
        if (!count || !create || !g.dgemm || !g.dmalloc || !g.dfree || !g.dmemcpy) return false;
        int n = 0;
        if (count(&n) != 0 || n < 1) return false;
        if (create(&g.handle) != 0) return false;
        g.devices = n;
        g.kind = "hip";
        g.name = "AMD GPU (hipBLAS, " + bl + ")";
        g.op_n = 111;  // HIPBLAS_OP_N
        g.op_t = 112;  // HIPBLAS_OP_T
        g.to_device = 1;  // hipMemcpyHostToDevice
        g.to_host = 2;    // hipMemcpyDeviceToHost
        return true;
    };
    if (mode == "cuda" || mode == "nvidia") g.ok = try_cuda();
    else if (mode == "hip" || mode == "amd" || mode == "rocm") g.ok = try_hip();
    else g.ok = try_cuda() || try_hip();
    return g;
}

struct State {
    CpuBlas blas;
    Gpu gpu;
    std::mutex gpu_mutex;
    std::atomic<std::size_t> gpu_products{0}, blas_products{0};
    std::atomic<double> gpu_flops{0.0};
};

State& state() {
    static State* s = [] {
        auto* st = new State;
        st->blas = load_cpu_blas();
        st->gpu = load_gpu();
        return st;
    }();
    return *s;
}

bool device_buffer(Gpu& g, int slot, std::size_t bytes) {
    if (g.cap[slot] >= bytes) return true;
    if (g.buf[slot]) g.dfree(g.buf[slot]);
    g.buf[slot] = nullptr;
    g.cap[slot] = 0;
    if (g.dmalloc(&g.buf[slot], bytes) != 0) return false;
    g.cap[slot] = bytes;
    return true;
}

// Elements spanned by a row-major matrix of `rows` rows with leading dimension ld and `cols`
// columns.
std::size_t extent(std::size_t rows, std::size_t cols, std::size_t ld) { return rows == 0 ? 0 : (rows - 1) * ld + cols; }

bool gpu_gemm(Gpu& g, bool ta, bool tb, std::size_t m, std::size_t n, std::size_t k, double alpha, const double* a,
              std::size_t lda, const double* b, std::size_t ldb, double beta, double* c, std::size_t ldc) {
    const std::size_t ea = ta ? extent(k, m, lda) : extent(m, k, lda);
    const std::size_t eb = tb ? extent(n, k, ldb) : extent(k, n, ldb);
    const std::size_t ec = extent(m, n, ldc);
    if (!device_buffer(g, 0, ea * sizeof(double)) || !device_buffer(g, 1, eb * sizeof(double)) ||
        !device_buffer(g, 2, ec * sizeof(double)))
        return false;
    if (g.dmemcpy(g.buf[0], a, ea * sizeof(double), g.to_device) != 0) return false;
    if (g.dmemcpy(g.buf[1], b, eb * sizeof(double), g.to_device) != 0) return false;
    if (beta != 0.0 && g.dmemcpy(g.buf[2], c, ec * sizeof(double), g.to_device) != 0) return false;
    // row-major C = op(A) op(B) is column-major C' = op(B)' op(A)': swap the operands
    const int st = g.dgemm(g.handle, tb ? g.op_t : g.op_n, ta ? g.op_t : g.op_n, static_cast<int>(n), static_cast<int>(m),
                           static_cast<int>(k), &alpha, static_cast<const double*>(g.buf[1]), static_cast<int>(ldb),
                           static_cast<const double*>(g.buf[0]), static_cast<int>(lda), &beta,
                           static_cast<double*>(g.buf[2]), static_cast<int>(ldc));
    if (st != 0) return false;
    return g.dmemcpy(c, g.buf[2], ec * sizeof(double), g.to_host) == 0;
}

}  // namespace

bool accelerated_gemm(bool trans_a, bool trans_b, std::size_t m, std::size_t n, std::size_t k, double alpha,
                      const double* a, std::size_t lda, const double* b, std::size_t ldb, double beta, double* c,
                      std::size_t ldc, int threads) {
    if (m == 0 || n == 0 || k == 0) return false;
    const std::size_t big = static_cast<std::size_t>(INT32_MAX);
    if (m > big || n > big || k > big || lda > big || ldb > big || ldc > big) return false;
    State& s = state();
    const double flops = 2.0 * static_cast<double>(m) * static_cast<double>(n) * static_cast<double>(k);
    if (s.gpu.ok && flops >= s.gpu.min_flops) {
        std::lock_guard<std::mutex> lock(s.gpu_mutex);
        if (gpu_gemm(s.gpu, trans_a, trans_b, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)) {
            ++s.gpu_products;
            double cur = s.gpu_flops.load();
            while (!s.gpu_flops.compare_exchange_weak(cur, cur + flops)) {
            }
            return true;
        }
        // a failed GPU product (e.g. out of device memory) falls through to the CPU
    }
    if (s.blas.g32 || s.blas.g64) {
        // small products stay with the portable code (no call overhead)
        if (flops < 2.0e5) return false;
        if (s.blas.set_threads) s.blas.set_threads(threads < 1 ? 1 : threads);
        if (s.blas.g32)
            s.blas.g32(kRowMajor, trans_a ? kTrans : kNoTrans, trans_b ? kTrans : kNoTrans, static_cast<int>(m),
                       static_cast<int>(n), static_cast<int>(k), alpha, a, static_cast<int>(lda), b, static_cast<int>(ldb),
                       beta, c, static_cast<int>(ldc));
        else
            s.blas.g64(kRowMajor, trans_a ? kTrans : kNoTrans, trans_b ? kTrans : kNoTrans, static_cast<std::int64_t>(m),
                       static_cast<std::int64_t>(n), static_cast<std::int64_t>(k), alpha, a, static_cast<std::int64_t>(lda),
                       b, static_cast<std::int64_t>(ldb), beta, c, static_cast<std::int64_t>(ldc));
        ++s.blas_products;
        return true;
    }
    return false;
}

std::string accelerator_description() {
    State& s = state();
    std::string d;
    if (s.gpu.ok) {
        char buf[64];
        std::snprintf(buf, sizeof buf, "%.0e", s.gpu.min_flops);
        d += s.gpu.name + ", " + std::to_string(s.gpu.devices) + " device" + (s.gpu.devices == 1 ? "" : "s") +
             ", for products of at least " + buf + " flop";
    }
    if (!s.blas.name.empty()) d += (d.empty() ? "" : "; ") + std::string("CPU BLAS ") + s.blas.name;
    return d;
}

AcceleratorCounters accelerator_counters() {
    State& s = state();
    AcceleratorCounters c;
    c.gpu_products = s.gpu_products.load();
    c.blas_products = s.blas_products.load();
    c.gpu_flops = s.gpu_flops.load();
    return c;
}

}  // namespace macqc
