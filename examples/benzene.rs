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
    println!("\n--- Structural Analysis ---");
    
    // 1. Bond Lengths
    println!("Bond Lengths (Å):");
    for i in 0..6 {
        let d_cc = (system.atoms[i].position - system.atoms[(i + 1) % 6].position).length();
        let d_ch = (system.atoms[i].position - system.atoms[i + 6].position).length();
        println!("  C{}-C{}: {:.4} Å | C{}-H{}: {:.4} Å", i+1, (i+1)%6 + 1, d_cc, i+1, i+7, d_ch);
    }

    // 2. Bond Angles
    println!("\nBond Angles (deg):");
    for j in 0..6 {
        let i = (j + 5) % 6;
        let k = (j + 1) % 6;
        let p_j = system.atoms[j].position;
        let p_i = system.atoms[i].position;
        let p_k = system.atoms[k].position;
        let v_ji = p_i - p_j;
        let v_jk = p_k - p_j;
        let angle = v_ji.angle_between(v_jk).to_degrees();
        println!("  C{}-C{}-C{}: {:.2}°", i+1, j+1, k+1, angle);
    }

    // 3. Dihedral Angles (Planarity)
    println!("\nDihedral Angles (deg) - for planarity check:");
    for i in 0..6 {
        let j = (i + 1) % 6;
        let k = (i + 2) % 6;
        let l = (i + 3) % 6;
        let p1 = system.atoms[i].position;
        let p2 = system.atoms[j].position;
        let p3 = system.atoms[k].position;
        let p4 = system.atoms[l].position;
        let b1 = p2 - p1; let b2 = p3 - p2; let b3 = p4 - p3;
        let n1 = b1.cross(b2).normalize();
        let n2 = b2.cross(b3).normalize();
        let m1 = n1.cross(b2.normalize());
        let x = n1.dot(n2); let y = m1.dot(n2);
        let phi = y.atan2(x).to_degrees();
        println!("  C{}-C{}-C{}-C{}: {:.2}°", i+1, j+1, k+1, l+1, phi);
    }
}
