use uff_relax::{System, Atom, Bond, UnitCell, UffOptimizer};
use glam::DVec3;

#[test]
fn test_polyethylene_pbc() {
    let mut atoms = Vec::new();
    let mut bonds = Vec::new();
    
    // Create a simple chain of 10 carbons
    let n_carbons = 10;
    for i in 0..n_carbons {
        // Add a small y displacement to avoid perfectly linear symmetry
        let pos = DVec3::new(i as f64 * 1.2, (i % 2) as f64 * 0.1, 0.0); 
        let atom = Atom::new(6, pos);
        atoms.push(atom);
        
        if i > 0 {
            bonds.push(Bond { atom_indices: (i-1, i), order: 1.0 });
        }
    }

    // Small box to force interactions with periodic images
    let cell = UnitCell::new_orthorhombic(DVec3::new(12.0, 5.0, 5.0));
    
    // Measure initial minimum distance (excluding bonded)
    let initial_min_dist = get_min_non_bonded_dist(&atoms, &bonds, &cell);
    println!("Initial min non-bonded distance: {:.4}", initial_min_dist);

    let mut system = System::new(atoms, bonds, cell);
    let optimizer = UffOptimizer::new(300, 1e-2).with_verbose(true);
    optimizer.optimize(&mut system);

    let final_min_dist = get_min_non_bonded_dist(&system.atoms, &system.bonds, &system.cell);
    println!("Final min non-bonded distance: {:.4}", final_min_dist);

    // In UFF, carbon vdW distance is ~3.8A, so it should be pushed away from its periodic image or neighbors
    assert!(final_min_dist > initial_min_dist, "Atoms should be pushed apart");
    assert!(final_min_dist > 2.0, "Overlaps should be resolved");
}

fn get_min_non_bonded_dist(atoms: &[Atom], bonds: &[Bond], cell: &UnitCell) -> f64 {
    let mut min_dist = f64::MAX;
    for i in 0..atoms.len() {
        for j in i+1..atoms.len() {
            let is_bonded = bonds.iter().any(|b| 
                (b.atom_indices.0 == i && b.atom_indices.1 == j) || 
                (b.atom_indices.0 == j && b.atom_indices.1 == i)
            );
            if !is_bonded {
                let d = cell.distance_vector(atoms[i].position, atoms[j].position).length();
                min_dist = min_dist.min(d);
            }
        }
    }
    min_dist
}
