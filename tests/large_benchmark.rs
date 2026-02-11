use uff_relax::{System, Atom, UnitCell, UffOptimizer};
use glam::DVec3;
use std::time::Instant;

fn create_large_system(n_atoms: usize) -> System {
    let mut atoms = Vec::new();
    let spacing = 2.5; 
    let side = (n_atoms as f64).powf(1.0/3.0).ceil() as usize;
    
    let mut count = 0;
    for x in 0..side {
        for y in 0..side {
            for z in 0..side {
                if count >= n_atoms { break; }
                let pos = DVec3::new(x as f64 * spacing, y as f64 * spacing, z as f64 * spacing);
                atoms.push(Atom::new(6, pos));
                count += 1;
            }
            if count >= n_atoms { break; }
        }
        if count >= n_atoms { break; }
    }

    let cell = UnitCell::new_orthorhombic(DVec3::new(500.0, 500.0, 500.0));
    System::new(atoms, Vec::new(), cell)
}

#[test]
fn test_default_optimization_100k() {
    let n = 100_000;
    println!("\nCreating system with {} atoms...", n);
    let mut system = create_large_system(n);
    
    // Default settings: threads=Auto (will be parallel), cutoff=6.0, history=10
    let optimizer = UffOptimizer::new(1000, 1.0).with_verbose(true);
    
    println!("Starting Optimization with default settings...");
    let start = Instant::now();
    optimizer.optimize(&mut system);
    let duration = start.elapsed();
    
    println!("\nTotal time for 100k atoms: {:?}", duration);
}