# uff-relax

A high-performance molecular geometry relaxation engine based on the **Universal Force Field (UFF)**, implemented in pure Rust.

Designed for rapid structural optimization of everything from small organic molecules to large-scale unit cells containing hundreds of thousands of atoms.

## Features

- **Full Periodic Table Support:** Implements parameters for elements H (1) to Lr (103) based on the original UFF paper (*Rappe et al. 1992*).
- **Parallel Processing:** Powered by `rayon` for efficient multi-threaded force calculations.
- **Large-Scale Optimization:** Uses a high-performance `CellList` implementation to maintain $O(N)$ scaling, enabling relaxation of systems with 100,000+ atoms in seconds.
- **Robust Convergence (FIRE):** Employs the **Fast Inertial Relaxation Engine (FIRE)** algorithm for stable and fast convergence even from highly distorted initial structures.
- **Periodic Boundary Conditions (PBC):** Full support for orthorhombic and triclinic unit cells.
- **Advanced Convergence Logic:** Combines $F_{\text{max}}$, $F_{\text{RMS}}$, and energy-stalling detection with rolling history averages to handle numerical noise and local instabilities.

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
uff-relax = "1.0.0"
```

## Quick Start

```rust
use uff_relax::{System, Atom, Bond, UnitCell, UffOptimizer};
use glam::DVec3;

fn main() {
    // 1. Define atoms and positions
    let atoms = vec![
        Atom::new(6, DVec3::new(0.0, 0.0, 0.0)), // Carbon
        Atom::new(1, DVec3::new(1.0, 0.0, 0.0)), // Hydrogen
    ];

    // 2. Define bonds
    let bonds = vec![
        Bond { atom_indices: (0, 1), order: 1.0 },
    ];

    // 3. Create a system (automatically infers UFF atom types)
    let cell = UnitCell::new_none();
    let mut system = System::new(atoms, bonds, cell);

    // 4. Run the optimizer
    let optimizer = UffOptimizer::new(1000, 0.1)
        .with_num_threads(4)
        .with_verbose(true);
    
    optimizer.optimize(&mut system);

    println!("Structural relaxation complete.");
}
```

## Performance & Scaling

`uff-relax` is optimized for modern hardware:
- **Small Molecules:** Automatically uses sequential processing to avoid threading overhead.
- **Large Systems:** Leverages `CellList` and parallel `rayon` kernels.
- **Continuity:** Implements a 5th-order polynomial switching function for non-bonded interactions to ensure smooth force transitions and reliable convergence.

## How it Works

1. **Parameter Assignment:** Automatically infers hybridization and UFF atom types based on element and connectivity.
2. **Spatial Partitioning:** Builds a dynamic `CellList` to accelerate non-bonded (Lennard-Jones) calculations.
3. **Multi-Objective Convergence:** Monitors the rolling average of maximum and root-mean-square forces to ensure the entire system is physically relaxed.

## License

MIT or Apache-2.0

## Author

**Forblaze Project**  
Website: [https://forblaze-works.com/](https://forblaze-works.com/)