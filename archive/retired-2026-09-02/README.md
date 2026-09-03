# Retired on 2026-09-02 (review finding F9)

Not compiled, not installed (`.Rbuildignore` excludes `archive/`). Kept for the record.

| File | Why retired |
|---|---|
| `src/matern_cholesky.cpp`, `src/matern_copula_cholesky.cpp` | Cholesky variants: no production caller; `marginal_sd_cholesky` materialises a dense N x N inverse (484 MB at 7,776 stations). `matern_copula_cholesky.cpp:58` carried the F3 bug (`D.log()` = Eigen's matrix logarithm; the archived copy has the one-line fix `D.array().log()`, verified against the dense reference before archiving). |
| `src/matern_copula_circulant.cpp`, `src/matern_copula_folded.cpp`, `src/circulant_utils.{cpp,h}` | FFTW circulant / folded approximations: no production caller; density, sampler and `construct_circulant_precision` measured mutually inconsistent (review §F9); per-call FFTW plans and buffers without RAII. Retiring them removes the fftw3 system dependency. |
| `R/plots.R` | `plot_matern()`: used only by the retired simulation diagnostics; pulled in undeclared ggplot2/dplyr imports. |
| `bench/`, `benchmark_all.csv` | The 2024 benchmark of the retired variants. The production-path benchmark is `tools/bench_production_paths.R`. |

Evidence: `notes/stdmatern-review-2026-09-02.html`.
