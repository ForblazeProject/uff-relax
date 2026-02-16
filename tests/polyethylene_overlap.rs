use uff_relax::{System, Atom, Bond, UnitCell, UffOptimizer};
use glam::DVec3;
use std::fs::File;
use std::io::{BufRead, BufReader};

fn parse_mol2(path: &str) -> (Vec<Atom>, Vec<Bond>, UnitCell) {
    let file = File::open(path).expect("Failed to open mol2 file");
    let reader = BufReader::new(file);
    let mut atoms = Vec::new();
    let mut bonds = Vec::new();
    let mut cell = UnitCell::new_none();
    
    let mut section = "";
    for line in reader.lines() {
        let line = line.unwrap();
        if line.starts_with("@<TRIPOS>ATOM") {
            section = "ATOM";
            continue;
        } else if line.starts_with("@<TRIPOS>BOND") {
            section = "BOND";
            continue;
        } else if line.starts_with("@<TRIPOS>CRYSIN") {
            section = "CRYSIN";
            continue;
        } else if line.starts_with("@<TRIPOS>") {
            section = "";
            continue;
        }

        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.is_empty() { continue; }

        match section {
            "ATOM" => {
                let element_str = parts[1].chars().next().unwrap();
                let element = match element_str {
                    'H' => 1,
                    'C' => 6,
                    _ => 6, // Default to Carbon for simplicity in this test
                };
                let x: f64 = parts[2].parse().unwrap();
                let y: f64 = parts[3].parse().unwrap();
                let z: f64 = parts[4].parse().unwrap();
                atoms.push(Atom::new(element, DVec3::new(x, y, z)));
            }
            "BOND" => {
                let u: usize = parts[1].parse::<usize>().unwrap() - 1;
                let v: usize = parts[2].parse::<usize>().unwrap() - 1;
                let order_str = parts[3];
                let order = match order_str {
                    "ar" => 1.5,
                    "2" => 2.0,
                    "3" => 3.0,
                    _ => 1.0,
                };
                bonds.push(Bond { atom_indices: (u, v), order });
            }
            "CRYSIN" => {
                let a: f64 = parts[0].parse().unwrap();
                let b: f64 = parts[1].parse().unwrap();
                let c: f64 = parts[2].parse().unwrap();
                cell = UnitCell::new_orthorhombic(DVec3::new(a, b, c));
            }
            _ => {}
        }
    }
    (atoms, bonds, cell)
}

#[test]
fn test_polyethylene_overlap_reproduction() {
    let (atoms, bonds, cell) = parse_mol2("pe_dist_cell_before_opt.mol2");
    let mut system = System::new(atoms, bonds, cell);

    for (i, atom) in system.atoms.iter().enumerate() {
        if uff_relax::get_uff_params(&atom.uff_type).is_none() {
            println!("Atom {} has type {} which has no UFF parameters!", i+1, atom.uff_type.as_str());
        }
    }

    let get_bad_pairs = |sys: &System, threshold: f64| {
        let n = sys.atoms.len();
        let mut adj = vec![Vec::new(); n];
        for b in &sys.bonds {
            adj[b.atom_indices.0].push(b.atom_indices.1);
            adj[b.atom_indices.1].push(b.atom_indices.0);
        }

        let mut bad_pairs = Vec::new();
        for i in 0..n {
            for j in i + 1..n {
                let mut is_bonded = false;
                for &n1 in &adj[i] {
                    if n1 == j { is_bonded = true; break; }
                    for &n2 in &adj[n1] {
                        if n2 == j { is_bonded = true; break; }
                    }
                }
                if !is_bonded {
                    let d = sys.cell.distance_vector(sys.atoms[i].position, sys.atoms[j].position).length();
                    if d < threshold { 
                        bad_pairs.push((i, j, d));
                    }
                }
            }
        }
        bad_pairs
    };

    let initial_bad = get_bad_pairs(&system, 1.0);
    println!("Initial bad pairs (< 1.0 A): {}", initial_bad.len());
    for (i, j, d) in initial_bad.iter().take(10) {
        println!("  Atom {} - {} : {:.6} A", i+1, j+1, d);
    }

    let optimizer = UffOptimizer::new(5000, 1e-3).with_verbose(true);
    optimizer.optimize(&mut system);

    let final_bad = get_bad_pairs(&system, 1.0);
    println!("Final bad pairs (< 1.0 A): {}", final_bad.len());
    for (i, j, d) in final_bad.iter().take(20) {
        println!("  Atom {} - {} : {:.6} A", i+1, j+1, d);
    }

    if !final_bad.is_empty() {
        let (i, j, d) = final_bad[0];
        panic!("Abnormal overlap remains: Atom {} - {} : {:.6} A", i+1, j+1, d);
    }
}
