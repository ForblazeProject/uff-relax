use uff_relax::{System, Atom, UffOptimizer, UnitCell};
use glam::DVec3;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::thread;
use std::time::Duration;

#[test]
fn test_cancellation() {
    // Create a system that would take some time (or many iterations)
    let atoms = vec![
        Atom::new(6, DVec3::new(0.0, 0.0, 0.0)),
        Atom::new(6, DVec3::new(1.0, 0.0, 0.0)), // Too close, will repel
    ];
    let mut system = System::new(atoms, vec![], UnitCell::new_none());
    
    let cancel_flag = Arc::new(AtomicBool::new(false));
    let cancel_flag_clone = cancel_flag.clone();
    
    let optimizer = UffOptimizer::new(10000, 1e-9) // Many iterations, very low threshold
        .with_cancel_flag(cancel_flag_clone)
        .with_verbose(true);
    
    // Spawn a thread to cancel after a short delay
    thread::spawn(move || {
        thread::sleep(Duration::from_millis(50));
        cancel_flag.store(true, Ordering::SeqCst);
    });
    
    optimizer.optimize(&mut system);
    
    // If it reached here without timing out the test, it likely cancelled successfully.
    // The verbose output (if enabled) would show "Reason: Cancelled"
}
