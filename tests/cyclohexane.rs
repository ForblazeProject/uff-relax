use uff_relax::{System, Atom, Bond, UnitCell, UffOptimizer};
use glam::DVec3;

#[test]
fn test_cyclohexane_chair() {
    let mut atoms = Vec::new();
    // 6 carbons in a slightly puckered ring to encourage chair
    for i in 0..6 {
        let angle = (i as f64) * 60.0f64.to_radians();
        let z = if i % 2 == 0 { 0.5 } else { -0.5 };
        atoms.push(Atom::new(6, DVec3::new(angle.cos() * 1.5, angle.sin() * 1.5, z)));
    }
    // 12 hydrogens (simplified, 2 per carbon)
    for i in 0..6 {
        atoms.push(Atom::new(1, atoms[i].position * 1.5));
        atoms.push(Atom::new(1, atoms[i].position * 1.5 + DVec3::new(0.0, 0.0, 0.5)));
    }

    let mut bonds = Vec::new();
    for i in 0..6 {
        bonds.push(Bond { atom_indices: (i, (i + 1) % 6), order: 1.0 });
        bonds.push(Bond { atom_indices: (i, 6 + i*2), order: 1.0 });
        bonds.push(Bond { atom_indices: (i, 6 + i*2 + 1), order: 1.0 });
    }

    let mut system = System::new(atoms, bonds, UnitCell::new_none());
    let optimizer = UffOptimizer::new(1000, 1e-3).with_verbose(true);
    optimizer.optimize(&mut system);

    // In a chair cyclohexane, dihedral angles are around +/- 50-60 deg
    println!("--- Cyclohexane Dihedrals ---");
    for i in 0..6 {
        let j = (i + 1) % 6; let k = (i + 2) % 6; let l = (i + 3) % 6;
        let phi = get_dihedral(&system, i, j, k, l);
        println!("Dihedral C{}-C{}-C{}-C{}: {:.2} deg", i, j, k, l, phi);
        assert!(phi.abs() > 40.0 && phi.abs() < 70.0);
    }
}

fn get_dihedral(system: &System, i: usize, j: usize, k: usize, l: usize) -> f64 {
    let p1 = system.atoms[i].position;
    let p2 = system.atoms[j].position;
    let p3 = system.atoms[k].position;
    let p4 = system.atoms[l].position;
    let b1 = p2 - p1; let b2 = p3 - p2; let b3 = p4 - p3;
    let n1 = b1.cross(b2).normalize();
    let n2 = b2.cross(b3).normalize();
    let m1 = n1.cross(b2.normalize());
    let x = n1.dot(n2); let y = m1.dot(n2);
    y.atan2(x).to_degrees()
}
