use uff_relax::{System, Atom, UnitCell, UffOptimizer};
use glam::DVec3;
use rand::Rng;

#[test]
fn test_random_high_density_clash() {
    let mut rng = rand::thread_rng();
    let n_atoms = 100;
    let l = 15.0; 
    let mut atoms = Vec::new();

    for _ in 0..n_atoms {
        let pos = DVec3::new(
            rng.gen_range(0.0..l),
            rng.gen_range(0.0..l),
            rng.gen_range(0.0..l),
        );
        atoms.push(Atom::new(18, pos)); // Argon (Z=18)
    }

    let cell = UnitCell::new_orthorhombic(DVec3::new(l, l, l));
    let mut system = System::new(atoms, Vec::new(), cell);

    let get_min_dist = |sys: &System| {
        let mut min_d = f64::MAX;
        for i in 0..sys.atoms.len() {
            for j in i + 1..sys.atoms.len() {
                let d = sys.cell.distance_vector(sys.atoms[i].position, sys.atoms[j].position).length();
                if d < min_d { min_d = d; }
            }
        }
        min_d
    };

    let initial_min = get_min_dist(&system);
    println!("Initial minimum distance: {:.6} A", initial_min);

    // Use a smaller step or more iterations to resolve heavy clashes
    let optimizer = UffOptimizer::new(2000, 1e-2).with_verbose(true);
    optimizer.optimize(&mut system);

    let final_min = get_min_dist(&system);
    println!("Final minimum distance: {:.6} A", final_min);

    // Ar vdW distance x1 is 3.868 A. min_dist should be around that.
    // We expect at least > 2.0 A for a stable relaxed structure.
    assert!(final_min > 2.0, "Abnormal overlap remains after relaxation: {:.6} A", final_min);
    assert!(final_min > initial_min, "Optimizer failed to increase minimum distance");
}
