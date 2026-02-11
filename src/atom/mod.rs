use glam::DVec3;
use serde::{Deserialize, Serialize};

/// Represents an atom type label in UFF (e.g., "C_3", "N_R").
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct UffAtomType(pub String);

impl UffAtomType {
    pub fn as_str(&self) -> &str {
        &self.0
    }
    pub fn unknown() -> Self {
        Self("Unknown".to_string())
    }
}

/// Represents a single atom in the system.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Atom {
    pub element: usize, // Atomic number
    pub position: DVec3,
    pub force: DVec3,
    pub uff_type: UffAtomType,
}

impl Atom {
    /// Creates a new atom with the given atomic number and position.
    pub fn new(element: usize, position: DVec3) -> Self {
        Self {
            element,
            position,
            force: DVec3::ZERO,
            uff_type: UffAtomType::unknown(),
        }
    }
}

/// Represents a chemical bond between two atoms.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Bond {
    pub atom_indices: (usize, usize),
    pub order: f32, // 1.0, 1.5, 2.0, 3.0
}
