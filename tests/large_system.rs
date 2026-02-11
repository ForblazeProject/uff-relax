use uff_relax::{System, Atom, UnitCell, UffOptimizer};
use glam::DVec3;

#[test]
fn test_large_system_parallel() {
    // Create a system with more than 500 atoms to trigger parallel path
    let mut atoms = Vec::new();
    let bonds = Vec::new();
    
    // Create 100 water-like molecules (300 atoms) - wait, need > 500
    // Let's create a 10x10x6 grid of Carbon atoms (600 atoms)
    let spacing = 2.0;
    for x in 0..10 {
        for y in 0..10 {
            for z in 0..6 {
                let pos = DVec3::new(x as f64 * spacing, y as f64 * spacing, z as f64 * spacing);
                atoms.push(Atom::new(6, pos));
            }
        }
    }

    let cell = UnitCell::new_orthorhombic(DVec3::new(20.0, 20.0, 12.0));
    let mut system = System::new(atoms, bonds, cell);
    
    let optimizer = UffOptimizer::new(5, 0.1).with_verbose(true);
    optimizer.optimize(&mut system);
    
    assert!(system.atoms.len() > 500);
}
