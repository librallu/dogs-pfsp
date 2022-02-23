//! DOGS implementation of the Graph Coloring problem

// #![warn(clippy::all, clippy::pedantic)]
// useful additional warnings if docs are missing, or crates imported but unused, etc.
#![warn(missing_debug_implementations)]
#![warn(missing_docs)]
#![warn(trivial_casts, trivial_numeric_casts)]
#![warn(unsafe_code)]
#![warn(unused_extern_crates)]
#![warn(variant_size_differences)]

// not sure if already by default in clippy
#![warn(clippy::similar_names)]
#![warn(clippy::shadow_unrelated)]
#![warn(clippy::shadow_same)]
#![warn(clippy::shadow_reuse)]

/// Defines permutation flowshop instances
pub mod pfsp;

/// forward flowtime minimization
pub mod fflowtime;

/// bidirectional makespan minimization
pub mod fbmakespan;

/// insertion-based branch-and-bound helper functions
pub mod nehhelper;

/// insertion-based branch-and-bound
pub mod nehmakespan;