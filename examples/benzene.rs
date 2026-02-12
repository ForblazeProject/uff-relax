use uff_relax::{System, Atom, Bond, UnitCell, UffOptimizer};
use glam::DVec3;

fn main() {
    // 1. Create atoms for Benzene (6 Carbons, 6 Hydrogens)
    let mut atoms = Vec::new();
    
    // Carbon ring (approximate positions)
    for i in 0..6 {
        let angle = (i as f64) * 60.0f64.to_radians();
        let pos = DVec3::new(angle.cos() * 1.4, angle.sin() * 1.4, 0.0);
        atoms.push(Atom::new(6, pos)); // 6 is Atomic Number for Carbon
    }
    
    // Hydrogen atoms (approximate positions)
    for i in 0..6 {
        let angle = (i as f64) * 60.0f64.to_radians();
        let pos = DVec3::new(angle.cos() * 2.4, angle.sin() * 2.4, 0.0);
        atoms.push(Atom::new(1, pos)); // 1 is Atomic Number for Hydrogen
    }

    // 2. Define bonds
    let mut bonds = Vec::new();
    for i in 0..6 {
        // C-C aromatic bonds (order 1.5)
        bonds.push(Bond { atom_indices: (i, (i + 1) % 6), order: 1.5 });
        // C-H single bonds (order 1.0)
        bonds.push(Bond { atom_indices: (i, i + 6), order: 1.0 });
    }

    // 3. Setup the system (No periodic boundary conditions)
    let cell = UnitCell::new_none();
    let mut system = System::new(atoms, bonds, cell);

    // 4. Configure the optimizer
    // Max 1000 iterations, force threshold 1e-3 kcal/mol/A
    let optimizer = UffOptimizer::new(1000, 1e-3)
        .with_verbose(true); // Print progress to console

    println!("Starting optimization of Benzene...");
    
    // 5. Run the optimization
    optimizer.optimize(&mut system);

    println!("Optimization complete!");
    
    // Check final C-C bond length
    let d01 = (system.atoms[0].position - system.atoms[1].position).length();
    println!("Final C-C bond length: {:.4} Å", d01);
}
