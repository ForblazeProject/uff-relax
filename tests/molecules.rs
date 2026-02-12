use uff_relax::{System, Atom, Bond, UnitCell, UffOptimizer};
use glam::DVec3;

fn get_dihedral_value(system: &System, i: usize, j: usize, k: usize, l: usize) -> f64 {
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

#[test]
fn test_benzene_optimization() {
    let mut atoms = Vec::new();
    // 6 carbons in a circle with some noise
    for i in 0..6 {
        let angle = (i as f64) * 60.0f64.to_radians();
        let pos = DVec3::new(angle.cos() * 1.3, angle.sin() * 1.3, 0.1 * (i % 2) as f64);
        let atom = Atom::new(6, pos);
        atoms.push(atom);
    }
    // 6 hydrogens
    for i in 0..6 {
        let angle = (i as f64) * 60.0f64.to_radians();
        let pos = DVec3::new(angle.cos() * 2.3, angle.sin() * 2.3, 0.0);
        let atom = Atom::new(1, pos);
        atoms.push(atom);
    }

    let mut bonds = Vec::new();
    for i in 0..6 {
        bonds.push(Bond { atom_indices: (i, (i + 1) % 6), order: 1.5 });
        bonds.push(Bond { atom_indices: (i, i + 6), order: 1.0 });
    }

    let cell = UnitCell::new_none();
    let mut system = System::new(atoms, bonds, cell);

    let optimizer = UffOptimizer::new(1000, 1e-3).with_verbose(true);
    optimizer.optimize(&mut system);

    // Check C-C bond length
    let d01 = (system.atoms[0].position - system.atoms[1].position).length();
    assert!((d01 - 1.40).abs() < 0.2);

    // Check CCC angles
    for j in 0..6 {
        let i = (j + 5) % 6;
        let k = (j + 1) % 6;
        let p_j = system.atoms[j].position;
        let p_i = system.atoms[i].position;
        let p_k = system.atoms[k].position;
        let v_ji = p_i - p_j;
        let v_jk = p_k - p_j;
        let cos_theta = v_ji.dot(v_jk) / (v_ji.length() * v_jk.length());
        let theta = cos_theta.acos().to_degrees();
        assert!((theta - 120.0).abs() < 15.0);
    }

    // Check planarity
    for i in 0..6 {
        let j = (i + 1) % 6; let k = (i + 2) % 6; let l = (i + 3) % 6;
        let phi = get_dihedral_value(&system, i, j, k, l);
        assert!(phi.abs() < 5.0);
    }
}

#[test]
fn test_cyclohexane_chair() {
    let mut atoms = Vec::new();
    // 6 carbons in a slightly puckered ring
    for i in 0..6 {
        let angle = (i as f64) * 60.0f64.to_radians();
        let z = if i % 2 == 0 { 0.5 } else { -0.5 };
        atoms.push(Atom::new(6, DVec3::new(angle.cos() * 1.5, angle.sin() * 1.5, z)));
    }
    // 12 hydrogens
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

    // Dihedral angles should be around +/- 50-60 deg
    for i in 0..6 {
        let j = (i + 1) % 6; let k = (i + 2) % 6; let l = (i + 3) % 6;
        let phi = get_dihedral_value(&system, i, j, k, l);
        assert!(phi.abs() > 40.0 && phi.abs() < 70.0);
    }
}

#[test]
fn test_polyethylene_pbc() {
    let mut atoms = Vec::new();
    let mut bonds = Vec::new();
    
    let n_carbons = 10;
    for i in 0..n_carbons {
        let pos = DVec3::new(i as f64 * 1.5, (i % 2) as f64 * 0.1, 0.0); 
        let atom = Atom::new(6, pos);
        atoms.push(atom);
        if i > 0 {
            bonds.push(Bond { atom_indices: (i-1, i), order: 1.0 });
        }
    }

    let cell = UnitCell::new_orthorhombic(DVec3::new(20.0, 10.0, 10.0));
    let mut system = System::new(atoms, bonds, cell);
    let optimizer = UffOptimizer::new(500, 1e-2).with_verbose(false);
    optimizer.optimize(&mut system);

    let final_min_dist = get_min_non_bonded_dist(&system.atoms, &system.bonds, &system.cell);
    // 1-4 interactions are scaled by 0.5, and C vdW is 3.85, so we expect some distance.
    // 1-3 are excluded. 1-4 should be around 2.5-3.0 A.
    assert!(final_min_dist > 2.0, "Non-bonded distance too small: {:.4}", final_min_dist);
}
