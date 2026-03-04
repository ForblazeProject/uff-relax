use uff_relax::{System, Atom, UffOptimizer, UnitCell};
use glam::DVec3;

#[test]
fn test_electrostatic_attraction() {
    // Two opposite charges
    let atoms = vec![
        Atom::new(6, DVec3::new(0.0, 0.0, 0.0)).with_charge(1.0),
        Atom::new(6, DVec3::new(5.0, 0.0, 0.0)).with_charge(-1.0),
    ];
    
    // Manual UFF types to avoid auto-assignment override if we want (actually auto-assign is fine)
    let mut system = System::new(atoms, vec![], UnitCell::new_none());
    
    // Initial energy check
    let energy_init = system.compute_forces_with_threads(1, 10.0);
    assert!(energy_init.electrostatic < 0.0, "Opposite charges should attract");
    
    // Optimize
    UffOptimizer::new(100, 1e-3)
        .with_cutoff(12.0)
        .optimize(&mut system);
    
    let dist = system.atoms[0].position.distance(system.atoms[1].position);
    // They should move closer together (though LJ repulsion will eventually stop them)
    assert!(dist < 5.0, "Charges should have moved closer. Final dist: {}", dist);
    
    let energy_final = system.compute_forces_with_threads(1, 10.0);
    assert!(energy_final.total < energy_init.total, "Energy should have decreased");
}

#[test]
fn test_electrostatic_repulsion() {
    // Two same charges
    let atoms = vec![
        Atom::new(6, DVec3::new(0.0, 0.0, 0.0)).with_charge(1.0),
        Atom::new(6, DVec3::new(2.0, 0.0, 0.0)).with_charge(1.0),
    ];
    
    let mut system = System::new(atoms, vec![], UnitCell::new_none());
    let energy_init = system.compute_forces_with_threads(1, 10.0);
    assert!(energy_init.electrostatic > 0.0, "Same charges should repel");
    
    // Optimize
    UffOptimizer::new(100, 1e-3)
        .with_cutoff(12.0)
        .optimize(&mut system);
    
    let dist = system.atoms[0].position.distance(system.atoms[1].position);
    assert!(dist > 2.0, "Charges should have moved apart. Final dist: {}", dist);
}
