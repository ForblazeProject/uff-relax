pub mod interactions;
pub mod parallel;
pub mod sequential;

use crate::atom::{Atom, Bond, UffAtomType};
use crate::cell::UnitCell;
use crate::params::element_symbol;
use glam::DVec3;

const PARALLEL_THRESHOLD: usize = 1000;

#[derive(Debug, Default, Clone, Copy)]
pub struct EnergyTerms {
    pub bond: f64,
    pub angle: f64,
    pub torsion: f64,
    pub non_bonded: f64,
    pub total: f64,
}

/// Represents a molecular system consisting of atoms, bonds, and an optional unit cell.
pub struct System {
    /// List of atoms in the system.
    pub atoms: Vec<Atom>,
    /// List of chemical bonds.
    pub bonds: Vec<Bond>,
    /// Unit cell for periodic boundary conditions.
    pub cell: UnitCell,
}

impl System {
    /// Creates a new molecular system and automatically assigns UFF atom types.
    ///
    /// # Arguments
    /// * `atoms` - Initial atom positions and elements.
    /// * `bonds` - Connectivity and bond orders.
    /// * `cell` - Boundary conditions (use `UnitCell::new_none()` for gas phase).
    pub fn new(atoms: Vec<Atom>, bonds: Vec<Bond>, cell: UnitCell) -> Self {
        let mut system = Self { atoms, bonds, cell };
        system.auto_assign_uff_types();
        system
    }

    /// Automatically infers UFF atom types based on element, connectivity, and bond orders.
    pub fn auto_assign_uff_types(&mut self) {
        let n = self.atoms.len();
        let mut adj = vec![Vec::new(); n];
        for bond in &self.bonds {
            adj[bond.atom_indices.0].push(bond);
            adj[bond.atom_indices.1].push(bond);
        }

        for i in 0..n {
            let z = self.atoms[i].element;
            let symbol = element_symbol(z);
            let neighbors = &adj[i];
            let n_neighbors = neighbors.len();
            let has_order_1_5 = neighbors.iter().any(|b| (b.order - 1.5).abs() < 0.1);
            let has_order_2_0 = neighbors.iter().any(|b| (b.order - 2.0).abs() < 0.1);
            let bond_order_sum: f32 = neighbors.iter().map(|b| b.order).sum();

            let label = match z {
                1 => "H_".to_string(),
                6 => { // Carbon
                    match n_neighbors {
                        4 => "C_3".to_string(),
                        3 => if has_order_1_5 || has_order_2_0 { "C_R".to_string() } else { "C_2".to_string() },
                        2 => "C_1".to_string(),
                        _ => "C_3".to_string(),
                    }
                }
                7 => { // Nitrogen
                    match n_neighbors {
                        3 => if has_order_1_5 { "N_R".to_string() } else { "N_3".to_string() },
                        2 => "N_2".to_string(),
                        1 => "N_1".to_string(),
                        _ => "N_3".to_string(),
                    }
                }
                8 => { // Oxygen
                    if has_order_1_5 { "O_R".to_string() }
                    else if n_neighbors == 1 && has_order_2_0 { "O_2".to_string() }
                    else { "O_3".to_string() }
                }
                9 => "F_".to_string(),
                15 => { // Phosphorus
                    if n_neighbors >= 4 || bond_order_sum > 4.0 { "P_3+5".to_string() }
                    else { "P_3+3".to_string() }
                }
                16 => { // Sulfur
                    if has_order_1_5 { "S_R".to_string() }
                    else if has_order_2_0 && n_neighbors == 1 { "S_2".to_string() }
                    else if n_neighbors == 3 || (bond_order_sum > 3.0 && bond_order_sum < 5.0) { "S_3+4".to_string() }
                    else if n_neighbors >= 4 || bond_order_sum >= 5.0 { "S_3+6".to_string() }
                    else { "S_3+2".to_string() }
                }
                17 => "Cl".to_string(),
                35 => "Br".to_string(),
                53 => "I_".to_string(),
                _ => {
                    if n_neighbors == 0 { format!("{}_", symbol) } 
                    else { format!("{}_", symbol) } // Always use symbol_ for generic match
                }
            };
            self.atoms[i].uff_type = UffAtomType(label);
        }
    }

    /// Computes forces and total energy breakdown.
    pub fn compute_forces(&mut self) -> EnergyTerms {
        self.compute_forces_with_threads(0, 6.0) // Default auto, cutoff 6.0
    }

    pub fn compute_forces_with_threads(&mut self, num_threads: usize, cutoff: f64) -> EnergyTerms {
        if num_threads == 1 {
            return self.compute_forces_serial(cutoff);
        }

        let use_parallel = if num_threads > 1 {
            true
        } else {
            self.atoms.len() >= PARALLEL_THRESHOLD
        };

        if use_parallel {
            let threads = if num_threads > 0 { num_threads } else { 4 };
            let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
            pool.install(|| self.compute_forces_parallel(cutoff))
        } else {
            self.compute_forces_serial(cutoff)
        }
    }

    fn compute_forces_serial(&mut self, cutoff: f64) -> EnergyTerms {
        let mut energy = EnergyTerms::default();
        for atom in &mut self.atoms { atom.force = DVec3::ZERO; }
        
        let mut adj = vec![Vec::new(); self.atoms.len()];
        for b in &self.bonds {
            let (u, v) = b.atom_indices;
            adj[u].push(v);
            adj[v].push(u);
        }

        energy.bond = self.compute_bond_forces_sequential();
        energy.angle = self.compute_angle_forces_sequential();
        energy.torsion = self.compute_torsion_forces_sequential();
        energy.non_bonded = self.compute_non_bonded_forces_sequential_cell_list(&adj, cutoff);
        energy.total = energy.bond + energy.angle + energy.torsion + energy.non_bonded;
        
        energy
    }

    fn compute_forces_parallel(&mut self, cutoff: f64) -> EnergyTerms {
        let mut energy = EnergyTerms::default();
        for atom in &mut self.atoms { atom.force = DVec3::ZERO; }
        
        let mut adj = vec![Vec::new(); self.atoms.len()];
        for b in &self.bonds {
            let (u, v) = b.atom_indices;
            adj[u].push(v);
            adj[v].push(u);
        }

        energy.bond = self.compute_bond_forces_parallel();
        energy.angle = self.compute_angle_forces_parallel();
        energy.torsion = self.compute_torsion_forces_parallel();
        energy.non_bonded = self.compute_non_bonded_forces_parallel_cell_list(&adj, cutoff);
        energy.total = energy.bond + energy.angle + energy.torsion + energy.non_bonded;
        
        energy
    }

    pub(crate) fn get_cell_neighbors(&self, cl: &crate::spatial::CellList, pos: DVec3, _cutoff: f64) -> Vec<usize> {
        let mut neighbors = Vec::new();
        let rel = pos - cl.min_p;
        let ix = (rel.x / cl.cell_size.x) as i32;
        let iy = (rel.y / cl.cell_size.y) as i32;
        let iz = (rel.z / cl.cell_size.z) as i32;

        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let mut nx = ix + dx; let mut ny = iy + dy; let mut nz = iz + dz;
                    
                    // PBC wrap for cells
                    if nx < 0 { nx += cl.dx as i32; } else if nx >= cl.dx as i32 { nx -= cl.dx as i32; }
                    if ny < 0 { ny += cl.dy as i32; } else if ny >= cl.dy as i32 { ny -= cl.dy as i32; }
                    if nz < 0 { nz += cl.dz as i32; } else if nz >= cl.dz as i32 { nz -= cl.dz as i32; }

                    if nx >= 0 && nx < cl.dx as i32 && ny >= 0 && ny < cl.dy as i32 && nz >= 0 && nz < cl.dz as i32 {
                        let idx = (nx as usize * cl.dy * cl.dz) + (ny as usize * cl.dz) + nz as usize;
                        neighbors.extend(&cl.cells[idx]);
                    }
                }
            }
        }
        
        // Remove duplicates (e.g. if cell size is large or system is small)
        neighbors.sort_unstable();
        neighbors.dedup();
        neighbors
    }

        pub(crate) fn get_exclusion_scale(&self, i: usize, j: usize, adj: &[Vec<usize>]) -> (bool, f64) {

            if i == j { return (true, 0.0); }

            

            // 1-2 neighbors: Strictly excluded from LJ

            for &n1 in &adj[i] {

                if n1 == j { return (true, 0.0); }

            }

            

            // 1-3 neighbors: Small scale (0.1) to prevent collapse

            for &n1 in &adj[i] {

                for &n2 in &adj[n1] {

                    if n2 == j { return (false, 0.1); }

                }

            }

            

            // 1-4 neighbors: Standard scale 0.5

            for &n1 in &adj[i] {

                for &n2 in &adj[n1] {

                    for &n3 in &adj[n2] {

                        if n3 == j { return (false, 0.5); }

                    }

                }

            }

            

            (false, 1.0)

        }

    }

    