/*
 * MacQC C API.
 *
 *   macqc_run(input_json)  runs HF (+ FCI) and returns the result as a JSON string.
 *                          It never throws and never returns NULL except when memory
 *                          is exhausted. Errors are reported in the "error" field.
 *   macqc_free(p)          releases a string returned by macqc_run.
 *
 * Jobs (cancellation and progress of long calculations, thread-safe):
 *   macqc_job_create()          a new job handle;
 *   macqc_run_job(json, job)    like macqc_run; the job may be cancelled from another thread,
 *                               then the result has "error": "calculation cancelled";
 *   macqc_job_cancel(job)       requests cancellation (returns at once);
 *   macqc_job_progress(job, …)  current stage (0 idle, 1 integrals, 2 SCF, 3 integral
 *                               transformation, 4 Full CI, 5 done, 6 coupled-cluster
 *                               iterations, 7 (T) correction), iteration (for 7: percent
 *                               done), total energy (E_h) and residual (SCF gradient, CI
 *                               residual or largest amplitude change);
 *                               (stages 8-13 of ΨLab 5.0: 8 gradient, 9 Hessian, 10 CPHF,
 *                               11 Lambda equations, 12 geometry optimisation, 13 numerical
 *                               Hessian);
 *   macqc_job_output(job, offset, &total)
 *                               the output file written so far from byte `offset` on (UTF-8,
 *                               free with macqc_free); `total` receives its full length, so
 *                               a reader polls with offset = the total of the previous call;
 *   macqc_job_progress_text(job) one line describing the current step, e.g.
 *                               "Optimisation step 4 · E = -76.0266 E_h · max force 2.1e-3"
 *                               (free with macqc_free);
 *   macqc_job_free(job)         releases the handle (after macqc_run_job has returned).
 *
 * Input keys of ΨLab 5.0 (H2, LiH, H2O): "job" (energy | opt | freq | opt+freq | opt+freq+sp),
 * "spMethod", "spBasis" (opt+freq+sp; default: method and basis), "temperature" (K, 298.15),
 * "pressure" (Pa, 101325), "optConvergence" (normal | tight), "maxSteps" (50). The result adds
 * "job", "optimization", "frequencies", "thermochemistry", "singlePoint" and "output" (the
 * complete output file); see README.
 *
 * Molecular orbitals (ΨLab 7.0). The result has "orbitals": the atoms, the basis shells and the
 * orbital sets ("sets": RHF canonical orbitals, all occupied and up to 10 virtual ones; Full CI
 * natural orbitals; for HeH+ in STO-nG the two Szabo-Ostlund orbitals) with their AO
 * coefficients. Pass that object back as text:
 *   macqc_orbital_surface(orbitals, set, orbital, isovalue, max_points, &positive, &negative)
 *                               the isosurfaces psi = +isovalue and psi = -isovalue (bohr;
 *                               normals point out of the lobes, triangles counter-clockwise
 *                               seen from outside) on a grid of at most max_points points
 *                               (0: 400000); release the meshes with macqc_mesh_free;
 *   macqc_orbital_plane(orbitals, set, orbital, origin, u, v, nu, nv, values)
 *                               the orbital at origin + i u + j v (bohr), i < nu fastest, j < nv,
 *                               into values (nu * nv floats).
 * Both return NULL on success, else an error message (free with macqc_free).
 */
#ifndef MACQC_API_H
#define MACQC_API_H

#ifdef __cplusplus
extern "C" {
#endif

char* macqc_run(const char* input_json);
void macqc_free(char* p);

typedef struct macqc_job macqc_job;
macqc_job* macqc_job_create(void);
char* macqc_run_job(const char* input_json, macqc_job* job);
void macqc_job_cancel(macqc_job* job);
void macqc_job_progress(const macqc_job* job, int* stage, int* iteration, double* energy, double* residual);
char* macqc_job_output(const macqc_job* job, long long offset, long long* total);
char* macqc_job_progress_text(const macqc_job* job);
void macqc_job_free(macqc_job* job);

typedef struct macqc_mesh {
    int vertex_count;
    float* positions;   /* 3 per vertex (bohr) */
    float* normals;     /* 3 per vertex, unit */
    int triangle_count;
    int* indices;       /* 3 per triangle */
} macqc_mesh;

char* macqc_orbital_surface(const char* orbitals_json, int set, int orbital, double isovalue, int max_points,
                            macqc_mesh* positive, macqc_mesh* negative);
char* macqc_orbital_plane(const char* orbitals_json, int set, int orbital, const double* origin, const double* u,
                          const double* v, int nu, int nv, float* values);
void macqc_mesh_free(macqc_mesh* mesh);

#ifdef __cplusplus
}
#endif

#endif /* MACQC_API_H */
