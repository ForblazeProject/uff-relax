use uff_relax::{System, Atom, UnitCell};
use glam::DVec3;

fn create_system(n_atoms: usize) -> System {
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
fn test_sequential_parallel_consistency() {
    let n = 2000;
    let mut system_seq = create_system(n);
    let mut system_par = create_system(n);
    let cutoff = 6.0;

    // Calculate forces sequentially
    let energy_seq = system_seq.compute_forces_with_threads(1, cutoff);
    let forces_seq: Vec<DVec3> = system_seq.atoms.iter().map(|a| a.force).collect();

    // Calculate forces in parallel
    let energy_par = system_par.compute_forces_with_threads(4, cutoff);
    let forces_par: Vec<DVec3> = system_par.atoms.iter().map(|a| a.force).collect();

    // Compare Energy
    let energy_diff = (energy_seq.total - energy_par.total).abs();
    println!("Energy Seq: {:.10}", energy_seq.total);
    println!("Energy Par: {:.10}", energy_par.total);
    println!("Energy Diff: {:.10e}", energy_diff);
    
    assert!(energy_diff < 1e-6, "Energy mismatch between sequential and parallel");

    // Compare Forces
    let mut max_force_diff = 0.0;
    for (f_s, f_p) in forces_seq.iter().zip(forces_par.iter()) {
        let diff = (*f_s - *f_p).length();
        if diff > max_force_diff {
            max_force_diff = diff;
        }
    }
    println!("Max Force Diff: {:.10e}", max_force_diff);
    assert!(max_force_diff < 1e-6, "Force mismatch between sequential and parallel");
}
