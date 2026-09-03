#ifndef STDMATERN_RNG_SEED_H
#define STDMATERN_RNG_SEED_H

#include <Rcpp.h>
#include <cstdint>
#include <random>

// Deterministic seeding for the eigen samplers.
//
// The samplers previously seeded one std::mt19937 per OpenMP thread from
// std::random_device -- which set.seed() cannot reach -- and shared a single
// std::normal_distribution across those threads (a data race on its cached
// Box-Muller value). Draws were therefore irreproducible from R.
//
// One generator per EIGENPAIR is seeded from (base seed, pair index), so the
// random stream is a function of the seed alone. Since 2026-09-02 the samplers
// run single-threaded (no OpenMP anywhere in the package), so the accumulation
// order is fixed too and a draw is bit-for-bit reproducible from the seed; the
// 2026-08-24 version still merged per-thread partial sums in arrival order and
// was reproducible only up to floating-point rounding at more than one thread.
// The seed space is 32 bits (std::mt19937 takes one word).

namespace stdmatern {

inline std::uint64_t splitmix64(std::uint64_t x) {
  x += 0x9E3779B97F4A7C15ULL;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31);
}

// The caller's seed when supplied; otherwise drawn from R's RNG so that
// set.seed() governs the draw (the Rcpp wrapper already opens an RNGScope).
inline std::uint64_t seed_base(Rcpp::Nullable<int> seed) {
  if (seed.isNotNull()) {
    return static_cast<std::uint64_t>(
        static_cast<std::uint32_t>(Rcpp::as<int>(seed)));
  }
  return static_cast<std::uint64_t>(
      static_cast<std::uint32_t>(R::unif_rand() * 4294967296.0));
}

// Generator for eigenpair (i, j) on a grid with dim_y columns.
inline std::mt19937 pair_generator(std::uint64_t base, int i, int j,
                                   int dim_y) {
  const std::uint64_t idx = static_cast<std::uint64_t>(i) * dim_y + j;
  return std::mt19937(
      static_cast<std::uint32_t>(splitmix64((base << 32) ^ idx)));
}

}  // namespace stdmatern

#endif  // STDMATERN_RNG_SEED_H
