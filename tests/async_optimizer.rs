use uff_relax::{System, Atom, UffOptimizer, UnitCell};
use glam::DVec3;
use futures::executor::block_on;
use std::sync::{Arc, Mutex};

#[test]
fn test_optimize_async() {
    let atoms = vec![
        Atom::new(6, DVec3::new(0.0, 0.0, 0.0)),
        Atom::new(6, DVec3::new(2.5, 0.0, 0.0)),
    ];
    let mut system = System::new(atoms, vec![], UnitCell::new_none());
    
    // Track if hook is called
    let hook_called = Arc::new(Mutex::new(0));
    let hook_called_clone = hook_called.clone();
    
    let optimizer = UffOptimizer::new(100, 1e-2)
        .with_step_hook(move |iter, _, _| {
            let mut count = hook_called_clone.lock().unwrap();
            *count += 1;
            if iter % 10 == 0 {
                // println!("Step: {}", iter);
            }
        });
    
    block_on(optimizer.optimize_async(&mut system));
    
    let dist = system.atoms[0].position.distance(system.atoms[1].position);
    // VdW for C is 3.85. 2.5 is too close, so they should move apart.
    assert!(dist > 2.5, "Atoms should have moved apart. Final dist: {}", dist);
    assert!(*hook_called.lock().unwrap() > 0, "Step hook should have been called");
}
