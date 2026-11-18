use std::collections::HashMap;
use std::f64::consts::PI;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::slice;

type V = [f64; 3];

#[derive(Clone)]
struct Rng(u64);

impl Rng {
    fn new(seed: u64) -> Self {
        Self(if seed == 0 { 0x9e37_79b9_7f4a_7c15 } else { seed })
    }

    fn next_u64(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / ((1u64 << 53) as f64)
    }

    fn index(&mut self, len: usize) -> usize {
        if len <= 1 {
            0
        } else {
            ((self.unit() * len as f64) as usize).min(len - 1)
        }
    }
}

fn add(a: V, b: V) -> V {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn sub(a: V, b: V) -> V {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn mul(a: V, x: f64) -> V {
    [a[0] * x, a[1] * x, a[2] * x]
}

fn dot(a: V, b: V) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: V, b: V) -> V {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(a: V) -> f64 {
    dot(a, a).sqrt()
}

fn unit(a: V) -> Option<V> {
    let length = norm(a);
    if !length.is_finite() || length < 1e-14 {
        None
    } else {
        Some(mul(a, 1.0 / length))
    }
}

fn rotate_point(point: V, axis: V, angle: f64, center: V) -> V {
    let Some(axis) = unit(axis) else {
        return point;
    };
    let v = sub(point, center);
    let (s, c) = angle.sin_cos();
    add(
        center,
        add(
            add(mul(v, c), mul(cross(axis, v), s)),
            mul(axis, dot(axis, v) * (1.0 - c)),
        ),
    )
}

fn perpendicular(v: V) -> V {
    let base = if v[0].abs() < 0.8 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
    unit(cross(v, base)).unwrap_or([0.0, 0.0, 1.0])
}

fn rotate_all(coords: &[V], axis: V, angle: f64, center: V) -> Vec<V> {
    coords.iter().map(|p| rotate_point(*p, axis, angle, center)).collect()
}

fn dihedral(p0: V, p1: V, p2: V, p3: V) -> Option<f64> {
    let b1 = sub(p1, p0);
    let b2 = sub(p2, p1);
    let b3 = sub(p3, p2);
    let n1 = cross(b1, b2);
    let n2 = cross(b2, b3);
    let b2n = unit(b2)?;
    let m1 = cross(n1, b2n);
    let x = dot(n1, n2);
    let y = dot(m1, n2);
    let angle = y.atan2(x);
    angle.is_finite().then_some(angle)
}

fn transform_unit(
    old: &[V],
    first: &[V],
    incoming: &[V],
    old_tail: usize,
    old_tail_neighbor: usize,
    first_head: usize,
    head: usize,
    tail: usize,
    head_neighbor: usize,
    bond_length: f64,
    target_dihedral: f64,
    random_rotation: bool,
    rng: &mut Rng,
) -> Option<Vec<V>> {
    let old_vec = sub(old[old_tail_neighbor], old[old_tail]);
    let new_vec = sub(incoming[head_neighbor], incoming[head]);
    let old_unit = unit(old_vec)?;
    let new_unit = unit(new_vec)?;
    let angle = dot(old_unit, new_unit).clamp(-1.0, 1.0).acos();
    let center = incoming[head_neighbor];

    // Rotate the incoming linker vector onto the continuation direction.
    // For parallel vectors, use a proper 180-degree rotation to preserve
    // molecular handedness.
    let (aligned, axis, align_angle) = if angle.abs() < 1e-12 {
        (incoming.to_vec(), perpendicular(new_unit), PI)
    } else if (PI - angle).abs() < 1e-12 {
        (incoming.to_vec(), [0.0, 0.0, 1.0], 0.0)
    } else {
        (incoming.to_vec(), cross(old_unit, new_unit), PI - angle)
    };
    let mut placed = rotate_all(&aligned, axis, align_angle, center);
    let anchor = sub(old[old_tail_neighbor], mul(old_unit, bond_length));
    let shift = sub(anchor, placed[head_neighbor]);
    for point in &mut placed {
        *point = add(*point, shift);
    }

    let twist = if random_rotation {
        rng.unit() * 2.0 * PI - PI
    } else {
        dihedral(
            first[first_head],
            old[old_tail_neighbor],
            placed[head_neighbor],
            placed[tail],
        )?
    };
    let placed_center = placed[head_neighbor];
    let result = rotate_all(&placed, mul(old_vec, -1.0), twist - target_dihedral, placed_center);
    if result.iter().flatten().any(|v| !v.is_finite()) {
        return None;
    }
    // The constructed inter-unit bond is expected to have the requested length.
    let bond = norm(sub(old[old_tail_neighbor], result[head_neighbor]));
    (bond.is_finite() && (bond - bond_length).abs() < 1e-8).then_some(result)
}

fn add_edge(adjacency: &mut [Vec<usize>], a: usize, b: usize) {
    adjacency[a].push(b);
    adjacency[b].push(a);
}

fn reachable_within_three(adjacency: &[Vec<usize>]) -> Vec<Vec<usize>> {
    let mut result = Vec::with_capacity(adjacency.len());
    let mut seen = vec![0usize; adjacency.len()];
    for start in 0..adjacency.len() {
        let stamp = start + 1;
        seen[start] = stamp;
        let mut found = Vec::new();
        let mut queue = vec![(start, 0usize)];
        let mut cursor = 0;
        while cursor < queue.len() {
            let (atom, depth) = queue[cursor];
            cursor += 1;
            if depth == 3 {
                continue;
            }
            for &neighbor in &adjacency[atom] {
                if seen[neighbor] != stamp {
                    seen[neighbor] = stamp;
                    found.push(neighbor);
                    queue.push((neighbor, depth + 1));
                }
            }
        }
        found.sort_unstable();
        result.push(found);
    }
    result
}

type Cell = (i64, i64, i64);

fn cell_of(point: V, size: f64) -> Cell {
    (
        (point[0] / size).floor() as i64,
        (point[1] / size).floor() as i64,
        (point[2] / size).floor() as i64,
    )
}

fn build_cells(
    placed: &[Vec<V>],
    atom_offsets: &[usize],
    heavy_atoms: &[Vec<usize>],
    cell_size: f64,
) -> HashMap<Cell, Vec<usize>> {
    let mut cells: HashMap<Cell, Vec<usize>> = HashMap::new();
    for (unit, coords) in placed.iter().enumerate() {
        let offset = atom_offsets[unit];
        for &local in &heavy_atoms[unit] {
            let point = coords[local];
            let global = offset + local;
            cells.entry(cell_of(point, cell_size)).or_default().push(global);
        }
    }
    cells
}

fn unit_internal_clear(
    unit: usize,
    coords: &[V],
    atom_offsets: &[usize],
    heavy_atoms: &[Vec<usize>],
    nearby: &[Vec<usize>],
    cutoff_squared: f64,
) -> bool {
    let global_offset = atom_offsets[unit];
    let atoms = &heavy_atoms[unit];
    for i in 0..atoms.len() {
        let a = atoms[i];
        for &b in &atoms[i + 1..] {
            if nearby[global_offset + a].binary_search(&(global_offset + b)).is_ok() {
                continue;
            }
            if dot(sub(coords[a], coords[b]), sub(coords[a], coords[b])) <= cutoff_squared {
                return false;
            }
        }
    }
    true
}


fn unit_external_clear(
    unit: usize,
    candidate: &[V],
    cells: &HashMap<Cell, Vec<usize>>,
    placed: &[Vec<V>],
    atom_offsets: &[usize],
    atom_to_unit: &[usize],
    heavy_atoms: &[Vec<usize>],
    nearby: &[Vec<usize>],
    cutoff_squared: f64,
    cell_size: f64,
) -> bool {
    let global_offset = atom_offsets[unit];
    for &local in &heavy_atoms[unit] {
        let point = candidate[local];
        let global = global_offset + local;
        let (cx, cy, cz) = cell_of(point, cell_size);
        for dx in -1i64..=1 {
            for dy in -1i64..=1 {
                for dz in -1i64..=1 {
                    let key = (cx.saturating_add(dx), cy.saturating_add(dy), cz.saturating_add(dz));
                    let Some(others) = cells.get(&key) else {
                        continue;
                    };
                    for &other in others {
                        if nearby[global].binary_search(&other).is_ok() {
                            continue;
                        }
                        let old_unit = atom_to_unit[other];
                        let old_local = other - atom_offsets[old_unit];
                        let delta = sub(point, placed[old_unit][old_local]);
                        if dot(delta, delta) <= cutoff_squared {
                            return false;
                        }
                    }
                }
            }
        }
    }
    true
}

fn bond_limit(order: f64) -> Option<f64> {
    if order == 1.0 {
        Some(2.7)
    } else if order == 1.5 {
        Some(1.9)
    } else if order == 2.0 {
        Some(1.8)
    } else if order == 3.0 {
        Some(1.4)
    } else {
        None
    }
}

fn bonds_valid(
    unit: usize,
    coords: &[V],
    heads: &[usize],
    tails: &[usize],
    bond_offsets: &[usize],
    bond_begin: &[usize],
    bond_end: &[usize],
    bond_orders: &[f64],
    check_bond_lengths: bool,
    skip_head: bool,
    skip_tail: bool,
) -> bool {
    if !check_bond_lengths {
        return true;
    }
    for bond in bond_offsets[unit]..bond_offsets[unit + 1] {
        let a = bond_begin[bond];
        let b = bond_end[bond];
        if (skip_head && (a == heads[unit] || b == heads[unit]))
            || (skip_tail && (a == tails[unit] || b == tails[unit]))
        {
            continue;
        }
        if let Some(limit) = bond_limit(bond_orders[bond]) {
            let delta = sub(coords[a], coords[b]);
            if dot(delta, delta) > limit * limit {
                return false;
            }
        }
    }
    true
}

fn read_unit_coords(
    coords: &[f64],
    coord_offsets: &[usize],
    atom_counts: &[usize],
    unit: usize,
    conformer: usize,
) -> Vec<V> {
    let start = coord_offsets[unit] + conformer * atom_counts[unit] * 3;
    let end = start + atom_counts[unit] * 3;
    coords[start..end]
        .chunks_exact(3)
        .map(|p| [p[0], p[1], p[2]])
        .collect()
}

struct GrowthResult {
    status: i32,
    coordinates: Vec<Vec<V>>,
    attempts: usize,
    rollbacks: usize,
    failed_step: usize,
}

fn grow_impl(
    unit_count: usize,
    total_atoms: usize,
    total_bonds: usize,
    total_coords: usize,
    coords: &[f64],
    coord_offsets: &[usize],
    atom_counts: &[usize],
    conformer_counts: &[usize],
    heads: &[usize],
    tails: &[usize],
    head_neighbors: &[usize],
    tail_neighbors: &[usize],
    atomic_numbers: &[u8],
    bond_offsets: &[usize],
    bond_begin: &[usize],
    bond_end: &[usize],
    bond_orders: &[f64],
    bond_length: f64,
    target_dihedral: f64,
    random_rotation: bool,
    dist_min: f64,
    check_bond_lengths: bool,
    retry_step: usize,
    retry: usize,
    rollback: usize,
    restart_initial: bool,
    seed: u64,
) -> GrowthResult {
    let invalid = || GrowthResult {
        status: 2,
        coordinates: Vec::new(),
        attempts: 0,
        rollbacks: 0,
        failed_step: 0,
    };
    if unit_count == 0
        || total_atoms == 0
        || retry_step == 0
        || !bond_length.is_finite()
        || bond_length <= 0.0
        || !target_dihedral.is_finite()
        || !dist_min.is_finite()
        || dist_min < 0.0
        || atom_counts.len() != unit_count
        || conformer_counts.len() != unit_count
        || heads.len() != unit_count
        || tails.len() != unit_count
        || head_neighbors.len() != unit_count
        || tail_neighbors.len() != unit_count
        || coord_offsets.len() != unit_count + 1
        || bond_offsets.len() != unit_count + 1
        || atomic_numbers.len() != total_atoms
        || bond_begin.len() != total_bonds
        || bond_end.len() != total_bonds
        || bond_orders.len() != total_bonds
    {
        return invalid();
    }

    let mut atom_offsets = Vec::with_capacity(unit_count + 1);
    atom_offsets.push(0usize);
    for &count in atom_counts {
        let Some(next) = atom_offsets.last().copied().unwrap().checked_add(count) else {
            return invalid();
        };
        atom_offsets.push(next);
    }
    if atom_offsets[unit_count] != total_atoms || coord_offsets[0] != 0 || bond_offsets[0] != 0 {
        return invalid();
    }
    for unit in 0..unit_count {
        let count = atom_counts[unit];
        if count == 0
            || conformer_counts[unit] == 0
            || heads[unit] >= count
            || tails[unit] >= count
            || heads[unit] == tails[unit]
            || head_neighbors[unit] >= count
            || tail_neighbors[unit] >= count
            || head_neighbors[unit] == heads[unit]
            || tail_neighbors[unit] == tails[unit]
            || head_neighbors[unit] == tails[unit]
            || tail_neighbors[unit] == heads[unit]
            || coord_offsets[unit + 1] < coord_offsets[unit]
            || bond_offsets[unit + 1] < bond_offsets[unit]
        {
            return invalid();
        }
        let Some(expected_coords) = conformer_counts[unit]
            .checked_mul(count)
            .and_then(|x| x.checked_mul(3))
        else {
            return invalid();
        };
        if coord_offsets[unit + 1] - coord_offsets[unit] != expected_coords {
            return invalid();
        }
        if bond_offsets[unit + 1] > total_bonds {
            return invalid();
        }
        for bond in bond_offsets[unit]..bond_offsets[unit + 1] {
            if bond_begin[bond] >= count
                || bond_end[bond] >= count
                || bond_begin[bond] == bond_end[bond]
                || !bond_orders[bond].is_finite()
            {
                return invalid();
            }
        }
    }
    if coord_offsets[unit_count] != total_coords || bond_offsets[unit_count] != total_bonds {
        return invalid();
    }
    if coords.iter().any(|x| !x.is_finite()) {
        return invalid();
    }

    let mut adjacency = vec![Vec::<usize>::new(); total_atoms];
    for unit in 0..unit_count {
        let offset = atom_offsets[unit];
        for bond in bond_offsets[unit]..bond_offsets[unit + 1] {
            let a = bond_begin[bond];
            let b = bond_end[bond];
            // Interior linker atoms are removed during assembly; the first
            // head and final tail remain as chain termini.
            let removed_head = unit > 0 && (a == heads[unit] || b == heads[unit]);
            let removed_tail = unit + 1 < unit_count && (a == tails[unit] || b == tails[unit]);
            if removed_head || removed_tail {
                continue;
            }
            add_edge(&mut adjacency, offset + a, offset + b);
        }
    }
    for unit in 1..unit_count {
        add_edge(
            &mut adjacency,
            atom_offsets[unit - 1] + tail_neighbors[unit - 1],
            atom_offsets[unit] + head_neighbors[unit],
        );
    }

    let mut atom_to_unit = vec![0usize; total_atoms];
    let mut heavy_atoms = vec![Vec::<usize>::new(); unit_count];
    for unit in 0..unit_count {
        for local in 0..atom_counts[unit] {
            let global = atom_offsets[unit] + local;
            atom_to_unit[global] = unit;
            let linker_is_removed = (unit > 0 && local == heads[unit])
                || (unit + 1 < unit_count && local == tails[unit]);
            if atomic_numbers[global] != 1 && !linker_is_removed {
                heavy_atoms[unit].push(local);
            }
        }
    }
    let nearby = reachable_within_three(&adjacency);
    let cutoff_squared = dist_min * dist_min;
    let cell_size = dist_min.max(1e-6);
    let mut growth_rng = Rng::new(seed ^ 0xa076_1d64_78bd_642f);
    let mut conformer_rng = Rng::new(seed.rotate_left(29) ^ 0xe703_7ed1_a0b4_28db);

    let first_conformer = if restart_initial {
        conformer_rng.index(conformer_counts[0])
    } else {
        0
    };
    let first_coords = read_unit_coords(
        coords,
        coord_offsets,
        atom_counts,
        0,
        first_conformer,
    );
    let mut placed = vec![first_coords];
    let mut base_valid = unit_internal_clear(
        0,
        &placed[0],
        &atom_offsets,
        &heavy_atoms,
        &nearby,
        cutoff_squared,
    ) && bonds_valid(
        0,
        &placed[0],
        heads,
        tails,
        bond_offsets,
        bond_begin,
        bond_end,
        bond_orders,
        check_bond_lengths,
        false,
        unit_count > 1,
    );

    if unit_count == 1 {
        return GrowthResult {
            status: 0,
            coordinates: placed,
            attempts: 0,
            rollbacks: 0,
            failed_step: 0,
        };
    }

    let mut history = vec![placed.clone()];
    let history_limit = rollback.saturating_add(1);
    let mut history_start_step = 0usize;
    let mut retries_left = retry;
    let mut attempts = 0usize;
    let mut rollbacks = 0usize;
    let mut start_step = 0usize;
    let mut failed_step = 0usize;

    while start_step < unit_count - 1 {
        let mut step = start_step;
        let mut success_chain = true;
        while step < unit_count - 1 {
            let unit = step + 1;
            let chosen_conformer = conformer_rng.index(conformer_counts[unit]);
            let incoming = read_unit_coords(
                coords,
                coord_offsets,
                atom_counts,
                unit,
                chosen_conformer,
            );
            let cells = build_cells(&placed, &atom_offsets, &heavy_atoms, cell_size);
            let mut accepted = None;

            for _ in 0..retry_step {
                attempts += 1;
                let Some(candidate) = transform_unit(
                    placed.last().unwrap(),
                    placed.first().unwrap(),
                    &incoming,
                    tails[step],
                    tail_neighbors[step],
                    heads[0],
                    heads[unit],
                    tails[unit],
                    head_neighbors[unit],
                    bond_length,
                    target_dihedral,
                    random_rotation,
                    &mut growth_rng,
                ) else {
                    continue;
                };

                let connection_bond_valid =
                    !check_bond_lengths || bond_length <= 2.7;
                let candidate_valid = base_valid
                    && connection_bond_valid
                    && unit_internal_clear(
                        unit,
                        &candidate,
                        &atom_offsets,
                        &heavy_atoms,
                        &nearby,
                        cutoff_squared,
                    )
                    && unit_external_clear(
                        unit,
                        &candidate,
                        &cells,
                        &placed,
                        &atom_offsets,
                        &atom_to_unit,
                        &heavy_atoms,
                        &nearby,
                        cutoff_squared,
                        cell_size,
                    )
                    && bonds_valid(
                        unit,
                        &candidate,
                        heads,
                        tails,
                        bond_offsets,
                        bond_begin,
                        bond_end,
                        bond_orders,
                        check_bond_lengths,
                        true,
                        false,
                    );
                if candidate_valid {
                    accepted = Some(candidate);
                    break;
                }
            }

            if let Some(candidate) = accepted {
                placed.push(candidate);
                history.push(placed.clone());
                if history.len() > history_limit {
                    history.remove(0);
                    history_start_step += 1;
                }
                step += 1;
            } else {
                failed_step = step + 1;
                success_chain = false;
                break;
            }
        }

        if success_chain {
            return GrowthResult {
                status: 0,
                coordinates: placed,
                attempts,
                rollbacks,
                failed_step: 0,
            };
        }

        if retries_left == 0 {
            return GrowthResult {
                status: 1,
                coordinates: Vec::new(),
                attempts,
                rollbacks,
                failed_step,
            };
        }

        retries_left -= 1;
        rollbacks += 1;
        let rollback_steps = rollback.min(step.saturating_sub(history_start_step));
        start_step = step - rollback_steps;
        let history_index = start_step - history_start_step;
        if history_index >= history.len() {
            return invalid();
        }
        placed = history[history_index].clone();
        history.truncate(history_index + 1);

        if start_step == 0 {
            let first_conformer = conformer_rng.index(conformer_counts[0]);
            placed[0] = read_unit_coords(
                coords,
                coord_offsets,
                atom_counts,
                0,
                first_conformer,
            );
            history[0] = placed.clone();
            base_valid = unit_internal_clear(
                0,
                &placed[0],
                &atom_offsets,
                &heavy_atoms,
                &nearby,
                cutoff_squared,
            ) && bonds_valid(
                0,
                &placed[0],
                heads,
                tails,
                bond_offsets,
                bond_begin,
                bond_end,
                bond_orders,
                check_bond_lengths,
                false,
                unit_count > 1,
            );
        } else {
            base_valid = true;
        }
    }

    GrowthResult {
        status: 0,
        coordinates: placed,
        attempts,
        rollbacks,
        failed_step: 0,
    }
}

unsafe fn raw_slice<'a, T>(ptr: *const T, len: usize) -> Option<&'a [T]> {
    if len == 0 {
        Some(&[])
    } else if ptr.is_null() {
        None
    } else {
        Some(slice::from_raw_parts(ptr, len))
    }
}

#[no_mangle]
pub unsafe extern "C" fn polymerize(
    unit_count: usize,
    total_coords: usize,
    total_atoms: usize,
    total_bonds: usize,
    coords_ptr: *const f64,
    coord_offsets_ptr: *const usize,
    atom_counts_ptr: *const usize,
    conformer_counts_ptr: *const usize,
    heads_ptr: *const usize,
    tails_ptr: *const usize,
    head_neighbors_ptr: *const usize,
    tail_neighbors_ptr: *const usize,
    atomic_numbers_ptr: *const u8,
    bond_offsets_ptr: *const usize,
    bond_begin_ptr: *const usize,
    bond_end_ptr: *const usize,
    bond_orders_ptr: *const f64,
    bond_length: f64,
    target_dihedral: f64,
    random_rotation: u8,
    dist_min: f64,
    check_bond_lengths: u8,
    retry_step: usize,
    retry: usize,
    rollback: usize,
    restart_initial: u8,
    seed: u64,
    output_ptr: *mut f64,
    stats_ptr: *mut usize,
) -> i32 {
    if unit_count == 0 || total_atoms == 0 || output_ptr.is_null() || stats_ptr.is_null() {
        return 2;
    }
    let Some(output_len) = total_atoms.checked_mul(3) else {
        return 2;
    };

    let result = catch_unwind(AssertUnwindSafe(|| {
        let Some(coords) = raw_slice(coords_ptr, total_coords) else { return 2; };
        let Some(coord_offsets) = raw_slice(coord_offsets_ptr, unit_count + 1) else { return 2; };
        let Some(atom_counts) = raw_slice(atom_counts_ptr, unit_count) else { return 2; };
        let Some(conformer_counts) = raw_slice(conformer_counts_ptr, unit_count) else { return 2; };
        let Some(heads) = raw_slice(heads_ptr, unit_count) else { return 2; };
        let Some(tails) = raw_slice(tails_ptr, unit_count) else { return 2; };
        let Some(head_neighbors) = raw_slice(head_neighbors_ptr, unit_count) else { return 2; };
        let Some(tail_neighbors) = raw_slice(tail_neighbors_ptr, unit_count) else { return 2; };
        let Some(atomic_numbers) = raw_slice(atomic_numbers_ptr, total_atoms) else { return 2; };
        let Some(bond_offsets) = raw_slice(bond_offsets_ptr, unit_count + 1) else { return 2; };
        let Some(bond_begin) = raw_slice(bond_begin_ptr, total_bonds) else { return 2; };
        let Some(bond_end) = raw_slice(bond_end_ptr, total_bonds) else { return 2; };
        let Some(bond_orders) = raw_slice(bond_orders_ptr, total_bonds) else { return 2; };

        let result = grow_impl(
            unit_count,
            total_atoms,
            total_bonds,
            total_coords,
            coords,
            coord_offsets,
            atom_counts,
            conformer_counts,
            heads,
            tails,
            head_neighbors,
            tail_neighbors,
            atomic_numbers,
            bond_offsets,
            bond_begin,
            bond_end,
            bond_orders,
            bond_length,
            target_dihedral,
            random_rotation != 0,
            dist_min,
            check_bond_lengths != 0,
            retry_step,
            retry,
            rollback,
            restart_initial != 0,
            seed,
        );

        let stats = slice::from_raw_parts_mut(stats_ptr, 3);
        stats[0] = result.attempts;
        stats[1] = result.rollbacks;
        stats[2] = result.failed_step;
        if result.status == 0 {
            let output = slice::from_raw_parts_mut(output_ptr, output_len);
            let mut cursor = 0usize;
            for unit in &result.coordinates {
                for point in unit {
                    output[cursor] = point[0];
                    output[cursor + 1] = point[1];
                    output[cursor + 2] = point[2];
                    cursor += 3;
                }
            }
            if cursor != output_len {
                return 2;
            }
        }
        result.status
    }));
    result.unwrap_or(3)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rng_repeats_for_same_seed() {
        let mut first = Rng::new(42);
        let mut second = Rng::new(42);
        for _ in 0..20 {
            assert_eq!(first.next_u64(), second.next_u64());
        }
    }

    #[test]
    fn transform_keeps_requested_bond_length() {
        let old = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
        let first = old.clone();
        let incoming = vec![[2.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let placed = transform_unit(
            &old,
            &first,
            &incoming,
            0,
            1,
            0,
            0,
            1,
            1,
            1.5,
            PI,
            true,
            &mut Rng::new(7),
        ).unwrap();
        assert!((norm(sub(old[1], placed[1])) - 1.5).abs() < 1e-8);
        assert!(placed.iter().flatten().all(|x| x.is_finite()));
    }
}
