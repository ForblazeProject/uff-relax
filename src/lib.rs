pub mod atom;
pub mod cell;
pub mod forcefield;
pub mod math;
pub mod optimizer;
pub mod params;
pub mod spatial;

pub use atom::{Atom, Bond, UffAtomType};
pub use cell::{UnitCell, CellType};
pub use forcefield::{System, EnergyTerms};
pub use optimizer::UffOptimizer;
pub use params::{get_uff_params, element_symbol};

use std::sync::Once;
static START: Once = Once::new();

/// Initializes the Rayon thread pool. 
/// If `num_threads` is Some(n), it sets that specific number.
/// If `num_threads` is None, it checks `RAYON_NUM_THREADS` env var or defaults to 4.
pub fn init_parallelism(num_threads: Option<usize>) {
    let threads = match num_threads {
        Some(n) => n,
        None => std::env::var("RAYON_NUM_THREADS")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(4),
    };

    // We can't re-initialize the global thread pool, so we only do it once.
    // For specific thread counts in benchmarks, we'll use a local thread pool
    // if we need to change it frequently, but for now we'll handle it via 
    // the compute_forces parameters.
    START.call_once(|| {
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global();
    });
}
