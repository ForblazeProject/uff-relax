use uff_relax::{System, Atom, UnitCell};
use glam::DVec3;
use std::time::Instant;

fn create_system(n_atoms: usize) -> System {
    let mut atoms = Vec::new();
    let spacing = 2.0;
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

    let cell = UnitCell::new_orthorhombic(DVec3::new(200.0, 200.0, 200.0));
    System::new(atoms, Vec::new(), cell)
}

fn measure_avg_time(system: &mut System, threads: usize, iters: usize, cutoff: f64) -> f64 {
    let start = Instant::now();
    for _ in 0..iters {
        system.compute_forces_with_threads(threads, cutoff);
    }
    start.elapsed().as_secs_f64() * 1000.0 / iters as f64
}

#[test]
fn test_scaling_benchmark() {
    let cutoff = 6.0;
    println!("\nBenchmark with Cutoff = {}", cutoff);
    println!("{:>10} | {:>12} | {:>12} | {:>8}", 
             "Atoms", "Serial (ms)", "Para-2 (ms)", "S-2x");
    println!("{:-<50}", "");

    for n in (100..=3000).step_by(100) {
        let mut system = create_system(n);
        
        // Warm up
        system.compute_forces_with_threads(1, cutoff);
        system.compute_forces_with_threads(2, cutoff);

        let t1 = measure_avg_time(&mut system, 1, 10, cutoff);
        let t2 = measure_avg_time(&mut system, 2, 10, cutoff);

        println!("{:>10} | {:>12.2} | {:>12.2} | {:>7.2}x", 
                 n, t1, t2, t1/t2);
    }
}