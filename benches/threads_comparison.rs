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

fn main() {
    let cutoff = 6.0;
    let thread_counts = [1, 2, 4, 8, 12, 16, 24, 32];
    let atom_counts = [1000, 5000, 10000];

    println!("
Thread Scalability Benchmark (Cutoff = {})", cutoff);
    println!("{:>10} | {:>10} | {:>10} | {:>8}", 
             "Atoms", "Threads", "Time (ms)", "Speedup");
    println!("{:-<50}", "");

    for &n in &atom_counts {
        let mut system = create_system(n);
        
        let t_serial = measure_avg_time(&mut system, 1, 5, cutoff);
        println!("{:>10} | {:>10} | {:>10.2} | {:>7.2}x", 
                 n, 1, t_serial, 1.0);

        for &threads in &thread_counts[1..] {
            let t = measure_avg_time(&mut system, threads, 5, cutoff);
            println!("{:>10} | {:>10} | {:>10.2} | {:>7.2}x", 
                     "", threads, t, t_serial / t);
        }
        println!("{:-<50}", "");
    }
}
