#pragma once

static constexpr bool USE_READ_CORRECTION = true;
static constexpr bool USE_COLUMN_MERGING  = true;
static constexpr bool MAP_TO_MAJORITY_INSTEAD_OF_N = false;
static constexpr int SUBSAMPLE_ITERS = 1000;

// Performance optimization settings
static constexpr int MAX_CACHED_EPP_SIZE = 2048;
static constexpr bool HIGH_MEMORY_CARTESIAN_MAP = true;
static constexpr int MUTEX_BIN_SIZE = 4096;
static constexpr int GRAIN_SIZE_FACTOR = 4;
static constexpr int NUM_RANGE_BINS = 50;
static constexpr int NUM_RANGE_TREES = 25;
static constexpr double SCORE_EPSILON = 1e-9;

// Internal Parameters
static constexpr double READ_DIST_FACTOR_THRESHOLD = (double) 0.5 / 100;
static constexpr int MAX_PEAK_PEAK_MUTATION = 2;
static constexpr int FREYJA_PEAKS_LIMIT = 5000;
static constexpr int TOP_N = 10;
static constexpr int MAX_PEAKS = 300;
static constexpr int MAX_NEIGHBORS_WEPP = 50;
static constexpr int MAX_NEIGHBORS_FREYJA = 500;
static constexpr int MAX_NEIGHBOR_MUTATION = 2;
// Update this to switch between Freyja and WEPP
static constexpr int MAX_NEIGHBOR_ITERATIONS = 1;
static constexpr bool FULL_TREE = false;

// Update based on real data or simulated data
static constexpr bool SIMULATED_DATA = false;

// Update if you are debugging and only want to run Iterative Freyja
static constexpr bool FULL_RUN = true;