# uff-relax

[![Crates.io](https://img.shields.io/crates/v/uff-relax.svg)](https://crates.io/crates/uff-relax)
[![Docs.rs](https://docs.rs/uff-relax/badge.svg)](https://docs.rs/uff-relax)
[![License](https://img.shields.io/badge/license-MIT%2FApache--2.0-blue.svg)](LICENSE-MIT)

A high-performance, parallelized molecular structure optimizer for Rust, powered by the **Universal Force Field (UFF)** and the **FIRE** algorithm.

## Features

- 🚀 **High Performance**: Optimized force evaluations with Cell Lists for efficient neighbor searching.
- 🧵 **Parallel Processing**: Scalable multi-threading via Rayon, automatically enabled for large systems (>1000 atoms) and defaulting to 4 threads for optimal performance.
- 💠 **PBC Support**: Periodic boundary conditions for Orthorhombic and Triclinic systems.
- 🧪 **Smart Type Assignment**: Automatically infers UFF atom types from atomic numbers and connectivity.
- 🦀 **Pure Rust**: Fast, safe, and easy to integrate.

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
uff-relax = "1.0.5" # Use the latest version
glam = { version = "0.31", features = ["serde"] }
rayon = "1.11"
serde = { version = "1.0", features = ["derive"] }
```

## Quick Start

```rust
use uff_relax::{System, Atom, Bond, UnitCell, UffOptimizer};
use glam::DVec3;

fn main() {
    // 1. Define Atoms
    let atoms = vec![
        Atom::new(6, DVec3::new(0.0, 0.0, 0.0)),
        Atom::new(6, DVec3::new(1.3, 0.0, 0.0)),
    ];

    // 2. Define Bonds
    let bonds = vec![
        Bond { atom_indices: (0, 1), order: 1.0 },
    ];

    // 3. Create System
    let mut system = System::new(atoms, bonds, UnitCell::new_none());

    // 4. Run Optimizer
    let optimizer = UffOptimizer::new(1000, 1e-3).with_verbose(true);
    optimizer.optimize(&mut system);

    println!("Final Energy: {:.4} kcal/mol", system.compute_forces().total);
}
```

## Running Examples

Try the included examples to see the optimizer in action:

```bash
cargo run --example benzene
```

## Benchmarks

This crate includes specialized benchmarks to measure scaling performance and handle large-scale systems. These are standalone binaries (`harness = false`) to ensure minimal overhead.

### 1. Parallelization Threshold (1 vs 2 threads)
Measures the efficiency of parallelization by comparing 1 and 2 threads as the number of atoms increases. This benchmark helped establish the initial `PARALLEL_THRESHOLD`.
```bash
cargo bench --bench threshold
```

### 2. Large System Stress Test
Simulates a system with 100,000 atoms to verify stability and memory efficiency in large-scale optimizations.
```bash
cargo bench --bench large_system_bench
```

### 3. Thread Scalability Comparison
Compares the performance across various thread counts (1, 2, 4, 8, etc.) for different system sizes, identifying optimal parallelization strategies.
```bash
cargo bench --bench threads_comparison
```

## License

Licensed under either of:

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
- MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.

## Author

**Forblaze Project**  
Website: [https://forblaze-works.com/](https://forblaze-works.com/)
