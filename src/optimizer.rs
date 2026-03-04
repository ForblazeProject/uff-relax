use tracing;
use crate::forcefield::System;
use glam::DVec3;
use web_time::Instant;
use std::collections::VecDeque;

/// Callback for optimization steps: (iteration, f_max, energy)
pub type StepHook = dyn Fn(usize, f64, f64) + Send + Sync;

/// Optimizer for molecular structures using the FIRE (Fast Iterative Relaxation Engine) algorithm.
pub struct UffOptimizer {
    /// Maximum number of iterations to perform.
    pub max_iterations: usize,
    /// Threshold for the maximum force on any atom (kcal/mol/Å).
    pub force_threshold: f64,
    /// Whether to print optimization progress to tracing logs.
    pub verbose: bool,
    /// Number of threads to use. 0 means automatic based on system size.
    pub num_threads: usize,
    /// Cutoff distance for non-bonded interactions (Å).
    pub cutoff: f64,
    /// Number of steps to average for convergence criteria.
    pub history_size: usize,
    /// Maximum distance an atom can move in a single step (Å).
    pub max_displacement: f64,
    /// Optional hook called after each iteration.
    pub step_hook: Option<std::sync::Arc<StepHook>>,
    /// Optional flag to cancel optimization from another thread/context.
    pub cancel_flag: Option<std::sync::Arc<std::sync::atomic::AtomicBool>>,
}

impl UffOptimizer {
    /// Creates a new optimizer with default settings.
    pub fn new(max_iterations: usize, force_threshold: f64) -> Self {
        Self {
            max_iterations,
            force_threshold,
            verbose: false,
            num_threads: 0,
            cutoff: 6.0,
            history_size: 10,
            max_displacement: 0.2,
            step_hook: None,
            cancel_flag: None,
        }
    }

    pub fn with_max_displacement(mut self, max: f64) -> Self {
        self.max_displacement = max;
        self
    }

    pub fn with_num_threads(mut self, num_threads: usize) -> Self {
        self.num_threads = num_threads;
        self
    }

    pub fn with_cutoff(mut self, cutoff: f64) -> Self {
        self.cutoff = cutoff;
        self
    }

    pub fn with_history_size(mut self, size: usize) -> Self {
        self.history_size = size;
        self
    }

    pub fn with_verbose(mut self, verbose: bool) -> Self {
        self.verbose = verbose;
        self
    }

    pub fn with_step_hook<F>(mut self, hook: F) -> Self 
    where F: Fn(usize, f64, f64) + Send + Sync + 'static {
        self.step_hook = Some(std::sync::Arc::new(hook));
        self
    }

    pub fn with_cancel_flag(mut self, flag: std::sync::Arc<std::sync::atomic::AtomicBool>) -> Self {
        self.cancel_flag = Some(flag);
        self
    }

    /// Optimized structural geometry using the FIRE algorithm.
    pub fn optimize(&self, system: &mut System) {
        let n = system.atoms.len();
        if n == 0 { return; }
        
        // Initial wrap only if periodic boundary conditions exist
        if !matches!(system.cell.cell_type, crate::cell::CellType::None) {
            for atom in &mut system.atoms {
                atom.position = system.cell.wrap_vector(atom.position);
            }
        }

        let mut velocities = vec![DVec3::ZERO; n];
        let mut dt = 0.01;
        let dt_max = 0.05;
        let mut n_pos = 0;
        let mut alpha = 0.15;
        let alpha_start = 0.15;

        let mut fmax_history = VecDeque::with_capacity(self.history_size);
        let mut frms_history = VecDeque::with_capacity(self.history_size);
        let mut ediff_history = VecDeque::with_capacity(self.history_size);
        let mut last_energy: Option<f64> = None;
        
        let start_time = Instant::now();

        if self.verbose {
            self.print_header(system);
        }

        let mut final_iter = 0;
        let mut final_status = "Max-Iter";

        for iter in 0..self.max_iterations {
            final_iter = iter;

            // Check for cancellation
            if let Some(ref cancel) = self.cancel_flag {
                if cancel.load(std::sync::atomic::Ordering::SeqCst) {
                    final_status = "Cancelled";
                    break;
                }
            }

            #[cfg(target_arch = "wasm32")]
            let energy = system.compute_forces_with_threads(1, self.cutoff);
            #[cfg(not(target_arch = "wasm32"))]
            let energy = system.compute_forces_with_threads(self.num_threads, self.cutoff);
            
            let (f_max, f_rms) = self.calculate_force_metrics(system);

            // Update history
            if fmax_history.len() >= self.history_size { fmax_history.pop_front(); }
            fmax_history.push_back(f_max);
            if frms_history.len() >= self.history_size { frms_history.pop_front(); }
            frms_history.push_back(f_rms);
            if let Some(prev_e) = last_energy {
                if ediff_history.len() >= self.history_size { ediff_history.pop_front(); }
                ediff_history.push_back((energy.total - prev_e).abs() / n as f64);
            }
            last_energy = Some(energy.total);

            // Convergence Check
            let (converged, status) = self.check_convergence(f_max, f_rms, &fmax_history, &frms_history, &ediff_history);
            
            if self.verbose && (iter % 10 == 0 || converged) {
                if energy.total.abs() >= 1e10 {
                    tracing::info!("{:>6} | {:>14.4} | {:>14.4} | {:>16.4e} | {:<10}", iter, f_max, f_rms, energy.total, status);
                } else {
                    tracing::info!("{:>6} | {:>14.4} | {:>14.4} | {:>16.4} | {:<10}", iter, f_max, f_rms, energy.total, status);
                }
            }

            if let Some(ref hook) = self.step_hook {
                hook(iter, f_max, energy.total);
            }

            if converged {
                final_status = status;
                break;
            }

            self.fire_update(system, &mut velocities, &mut dt, dt_max, &mut n_pos, &mut alpha, alpha_start);
        }

        if self.verbose {
            self.print_footer(system, final_status, start_time, final_iter, &fmax_history, &frms_history);
        }
    }

    /// Asynchronous version of the optimizer for non-blocking environments (Wasm/UIs).
    pub async fn optimize_async(&self, system: &mut System) {
        let n = system.atoms.len();
        if n == 0 { return; }
        
        if !matches!(system.cell.cell_type, crate::cell::CellType::None) {
            for atom in &mut system.atoms {
                atom.position = system.cell.wrap_vector(atom.position);
            }
        }

        let mut velocities = vec![DVec3::ZERO; n];
        let mut dt = 0.01;
        let dt_max = 0.05;
        let mut n_pos = 0;
        let mut alpha = 0.15;
        let alpha_start = 0.15;

        let mut fmax_history = VecDeque::with_capacity(self.history_size);
        let mut frms_history = VecDeque::with_capacity(self.history_size);
        let mut ediff_history = VecDeque::with_capacity(self.history_size);
        let mut last_energy: Option<f64> = None;
        
        let start_time = Instant::now();

        if self.verbose {
            self.print_header(system);
        }

        let mut final_iter = 0;
        let mut final_status = "Max-Iter";

        for iter in 0..self.max_iterations {
            final_iter = iter;

            // Check for cancellation
            if let Some(ref cancel) = self.cancel_flag {
                if cancel.load(std::sync::atomic::Ordering::SeqCst) {
                    final_status = "Cancelled";
                    break;
                }
            }
            
            // In Wasm, num_threads should be 1 as Rayon is not supported easily
            #[cfg(target_arch = "wasm32")]
            let energy = system.compute_forces_with_threads(1, self.cutoff);
            #[cfg(not(target_arch = "wasm32"))]
            let energy = system.compute_forces_with_threads(self.num_threads, self.cutoff);
            
            let (f_max, f_rms) = self.calculate_force_metrics(system);

            if fmax_history.len() >= self.history_size { fmax_history.pop_front(); }
            fmax_history.push_back(f_max);
            if frms_history.len() >= self.history_size { frms_history.pop_front(); }
            frms_history.push_back(f_rms);
            if let Some(prev_e) = last_energy {
                if ediff_history.len() >= self.history_size { ediff_history.pop_front(); }
                ediff_history.push_back((energy.total - prev_e).abs() / n as f64);
            }
            last_energy = Some(energy.total);

            let (converged, status) = self.check_convergence(f_max, f_rms, &fmax_history, &frms_history, &ediff_history);
            
            if self.verbose && (iter % 10 == 0 || converged) {
                if energy.total.abs() >= 1e10 {
                    tracing::info!("{:>6} | {:>14.4} | {:>14.4} | {:>16.4e} | {:<10}", iter, f_max, f_rms, energy.total, status);
                } else {
                    tracing::info!("{:>6} | {:>14.4} | {:>14.4} | {:>16.4} | {:<10}", iter, f_max, f_rms, energy.total, status);
                }
            }

            if let Some(ref hook) = self.step_hook {
                hook(iter, f_max, energy.total);
            }

            if converged {
                final_status = status;
                break;
            }

            self.fire_update(system, &mut velocities, &mut dt, dt_max, &mut n_pos, &mut alpha, alpha_start);

            // Yield control back to the environment periodically
            if iter % 5 == 0 {
                self.yield_now().await;
            }
        }

        if self.verbose {
            self.print_footer(system, final_status, start_time, final_iter, &fmax_history, &frms_history);
        }
    }

    async fn yield_now(&self) {
        #[cfg(feature = "wasm")]
        {
            let promise = js_sys::Promise::new(&mut |resolve, _| {
                if let Some(window) = web_sys::window() {
                    window.set_timeout_with_callback_and_timeout_and_arguments_0(&resolve, 0).unwrap();
                }
            });
            let _ = wasm_bindgen_futures::JsFuture::from(promise).await;
        }
    }

    fn calculate_force_metrics(&self, system: &System) -> (f64, f64) {
        let n = system.atoms.len();
        let mut max_f_sq: f64 = 0.0;
        let mut sum_f_sq: f64 = 0.0;
        for atom in &system.atoms {
            let f_sq = atom.force.length_squared();
            max_f_sq = f64::max(max_f_sq, f_sq);
            sum_f_sq += f_sq;
        }
        (max_f_sq.sqrt(), (sum_f_sq / (3.0 * n as f64)).sqrt())
    }

    fn check_convergence(&self, _f_max: f64, _f_rms: f64, fmax_hist: &VecDeque<f64>, frms_hist: &VecDeque<f64>, ediff_hist: &VecDeque<f64>) -> (bool, &'static str) {
        if fmax_hist.len() < self.history_size {
            return (false, "");
        }
        let avg_fmax: f64 = fmax_hist.iter().sum::<f64>() / self.history_size as f64;
        let avg_frms: f64 = frms_hist.iter().sum::<f64>() / self.history_size as f64;
        let avg_ediff: f64 = if ediff_hist.is_empty() { 1.0 } else { ediff_hist.iter().sum::<f64>() / ediff_hist.len() as f64 };

        if avg_fmax < self.force_threshold {
            (true, "Fmax-Conv")
        } else if avg_fmax < self.force_threshold * 2.0 && avg_frms < self.force_threshold * 0.5 {
            (true, "FRMS-Conv")
        } else if !ediff_hist.is_empty() && avg_ediff < 1e-7 {
            (true, "E-Stalled")
        } else {
            (false, "")
        }
    }

    fn fire_update(&self, system: &mut System, velocities: &mut [DVec3], dt: &mut f64, dt_max: f64, n_pos: &mut usize, alpha: &mut f64, alpha_start: f64) {
        let n = system.atoms.len();
        let mut p = 0.0;
        for i in 0..n {
            p += velocities[i].dot(system.atoms[i].force);
        }

        for i in 0..n {
            let f_norm = system.atoms[i].force.length();
            let v_norm = velocities[i].length();
            if f_norm > 1e-9 {
                velocities[i] = (1.0 - *alpha) * velocities[i] + *alpha * (system.atoms[i].force / f_norm) * v_norm;
            }
        }

        if p > 0.0 {
            *n_pos += 1;
            if *n_pos > 5 {
                *dt = f64::min(*dt * 1.05, dt_max);
                *alpha *= 0.99;
            }
        } else {
            *n_pos = 0;
            *dt *= 0.5;
            *alpha = alpha_start;
            for v in velocities.iter_mut() {
                *v = DVec3::ZERO;
            }
        }

        for i in 0..n {
            velocities[i] += system.atoms[i].force * (*dt);
            let mut move_vec = velocities[i] * (*dt);
            let move_len = move_vec.length();
            if move_len > self.max_displacement {
                move_vec *= self.max_displacement / move_len;
                velocities[i] = move_vec / (*dt);
            }
            let new_pos = system.atoms[i].position + move_vec;
            system.atoms[i].position = system.cell.wrap_vector(new_pos);
        }
    }

    fn print_header(&self, system: &System) {
        let n_atoms = system.atoms.len();
        let n_bonds = system.bonds.len();
        let has_charges = system.atoms.iter().any(|a| a.charge.abs() > 1e-12);
        
        // Determine actual threads used
        #[cfg(target_arch = "wasm32")]
        let actual_threads = 1;
        
        #[cfg(not(target_arch = "wasm32"))]
        let actual_threads = if self.num_threads == 1 {
            1
        } else if self.num_threads > 1 {
            self.num_threads
        } else if n_atoms >= 1000 { // PARALLEL_THRESHOLD
            std::env::var("RAYON_NUM_THREADS")
                .ok()
                .and_then(|s| s.parse().ok())
                .unwrap_or(4)
        } else {
            1
        };

        let version_str = format!(" uff-relax v{} ", env!("CARGO_PKG_VERSION"));
        tracing::info!("\n{:=^80}", version_str);
        tracing::info!("{:<10} {:<10} | {:<10} {:<10}", "Atoms:", n_atoms, "Bonds:", n_bonds);
        tracing::info!("{:<10} {:<10.1} | {:<10} {:<10.4} kcal/mol/Å", "Cutoff:", self.cutoff, "Threshold:", self.force_threshold);
        tracing::info!("{:<10} {:<10} | {:<10} {:<10}", 
            "Threads:", actual_threads, 
            "Charges:", if has_charges { "Active (Wolf)" } else { "Inactive" }
        );
        tracing::info!("{:<10} {:<10} | {:<10} {:<10}", "Max Iter:", self.max_iterations, "", "");
        tracing::info!("{:-<80}", "");
        tracing::info!("{:<6} | {:<14} | {:<14} | {:<16} | {:<10}", "", "Fmax", "FRMS", "Total E", "");
        tracing::info!("{:<6} | {:<14} | {:<14} | {:<16} | {:<10}", "Iter", "(kcal/mol/Å)", "(kcal/mol/Å)", "(kcal/mol)", "Status");
        tracing::info!("{:-<80}", "");
    }

    fn print_footer(&self, system: &mut System, final_status: &str, start_time: Instant, final_iter: usize, fmax_hist: &VecDeque<f64>, frms_hist: &VecDeque<f64>) {
        let n = system.atoms.len();
        let duration = start_time.elapsed();
        let final_energy = system.compute_forces_with_threads(self.num_threads, self.cutoff);

        let mut min_dist = f64::MAX;
        let mut min_pair = (0, 0);
        for i in 0..n {
            for j in i + 1..n {
                let d = system.cell.distance_vector(system.atoms[i].position, system.atoms[j].position).length();
                if d < min_dist { 
                    min_dist = d;
                    min_pair = (i, j);
                }
            }
        }

        tracing::info!("{:-<80}", "");
        tracing::info!("=== Optimization Finished ===");
        tracing::info!("Reason: {:<20}", final_status);
        tracing::info!("Total Time: {:<10.3?} (Avg: {:.3?} / step)", duration, duration / (final_iter + 1) as u32);
        if final_energy.total.abs() >= 1e10 {
            tracing::info!("Final Energy: {:<15.4e} kcal/mol", final_energy.total);
        } else {
            tracing::info!("Final Energy: {:<15.4} kcal/mol", final_energy.total);
        }
        tracing::info!("Final Fmax:   {:<15.4} kcal/mol/Å", fmax_hist.back().unwrap_or(&0.0));
        tracing::info!("Final FRMS:   {:<15.4} kcal/mol/Å", frms_hist.back().unwrap_or(&0.0));
        tracing::info!("Min Distance: {:<15.4} Å (Atoms {} and {})", min_dist, min_pair.0 + 1, min_pair.1 + 1);
        tracing::info!("{:>80}", "(c) 2026 Forblaze Project");
        tracing::info!("{:-<80}\n", "");
    }
}