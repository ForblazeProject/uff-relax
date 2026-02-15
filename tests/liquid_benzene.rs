use uff_relax::{System, Atom, Bond, UnitCell, UffOptimizer};
use glam::DVec3;

fn create_benzene(offset: DVec3, start_idx: usize) -> (Vec<Atom>, Vec<Bond>) {
    let mut atoms = Vec::new();
    let mut bonds = Vec::new();
    
    // 6 carbons
    for i in 0..6 {
        let angle = (i as f64) * 60.0f64.to_radians();
        let pos = DVec3::new(angle.cos() * 1.4, angle.sin() * 1.4, 0.0) + offset;
        atoms.push(Atom::new(6, pos));
    }
    // 6 hydrogens
    for i in 0..6 {
        let angle = (i as f64) * 60.0f64.to_radians();
        let pos = DVec3::new(angle.cos() * 2.5, angle.sin() * 2.5, 0.0) + offset;
        atoms.push(Atom::new(1, pos));
    }
    // Bonds
    for i in 0..6 {
        bonds.push(Bond { atom_indices: (start_idx + i, start_idx + (i + 1) % 6), order: 1.5 });
        bonds.push(Bond { atom_indices: (start_idx + i, start_idx + i + 6), order: 1.0 });
    }
    (atoms, bonds)
}

#[test]
fn test_benzene_liquid_density() {
    // Benzene MW = 78.11 g/mol
    // Target density = 0.8 g/cm3
    // For 8 molecules:
    // Volume = (8 * 78.11) / (0.8 * 0.6022) = 1297 A3
    // Side length L = 10.9 A
    let _n_mols = 8;
    let l = 10.9;
    let mut all_atoms = Vec::new();
    let mut all_bonds = Vec::new();

    for i in 0..2 {
        for j in 0..2 {
            for k in 0..2 {
                let offset = DVec3::new(i as f64 * l/2.0 + l/4.0, j as f64 * l/2.0 + l/4.0, k as f64 * l/2.0 + l/4.0);
                let (atoms, bonds) = create_benzene(offset, all_atoms.len());
                all_atoms.extend(atoms);
                all_bonds.extend(bonds);
            }
        }
    }

    let cell = UnitCell::new_orthorhombic(DVec3::new(l, l, l));
    let mut system = System::new(all_atoms, all_bonds, cell);

    // Initial check: some atoms might be close due to simple grid placement
    let optimizer = UffOptimizer::new(1000, 1e-2).with_verbose(false);
    optimizer.optimize(&mut system);

    // Final check for abnormal distances
    let mut min_dist = 100.0f64;
    let n = system.atoms.len();
    for i in 0..n {
        for j in i + 1..n {
            let d = system.cell.distance_vector(system.atoms[i].position, system.atoms[j].position).length();
            if d < min_dist {
                min_dist = d;
            }
        }
    }

    println!("Minimum distance in relaxed liquid benzene: {:.4} A", min_dist);
    // Even for bonded H-C, distance is ~1.1 A. 
    // For non-bonded, it should be > 1.8 A (H-H vdW is 2.88, so 0.4*2.88 is 1.15)
    // If it's 0.75 A or 0.1 A, it's a bug.
    assert!(min_dist > 0.95, "Abnormal atomic overlap detected: {:.4} A", min_dist);
}
