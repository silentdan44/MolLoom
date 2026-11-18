"""ctypes adapter for the optional Rust polymer growth core."""

from __future__ import annotations

import ctypes
from functools import lru_cache
import os
from pathlib import Path
import sys

import numpy as np


def _library_name() -> str:
    if os.name == "nt":
        return "polymer_growth_core.dll"
    if sys.platform == "darwin":
        return "libpolymer_growth_core.dylib"
    return "libpolymer_growth_core.so"


@lru_cache(maxsize=1)
def _core():
    name = _library_name()
    package_dir = Path(__file__).resolve().parent
    candidates = [
        package_dir / name,
        package_dir / "rust_growth" / "target" / "release" / name,
        package_dir / "rust_growth" / "target" / "debug" / name,
    ]
    library_path = next((path for path in candidates if path.is_file()), None)
    if library_path is None:
        raise ImportError(
            "Rust growth core is not built. From the project root, run "
            "'cargo build --release --manifest-path polymer_lib/rust_growth/Cargo.toml'."
        )

    library = ctypes.CDLL(str(library_path))
    function = library.polymerize
    size_ptr = ctypes.POINTER(ctypes.c_size_t)
    double_ptr = ctypes.POINTER(ctypes.c_double)
    u8_ptr = ctypes.POINTER(ctypes.c_uint8)
    function.argtypes = [
        ctypes.c_size_t, ctypes.c_size_t, ctypes.c_size_t, ctypes.c_size_t,
        double_ptr, size_ptr, size_ptr, size_ptr, size_ptr, size_ptr,
        size_ptr, size_ptr, u8_ptr, size_ptr, size_ptr, size_ptr, double_ptr,
        ctypes.c_double, ctypes.c_double, ctypes.c_uint8, ctypes.c_double,
        ctypes.c_uint8, ctypes.c_size_t, ctypes.c_size_t, ctypes.c_size_t,
        ctypes.c_uint8, ctypes.c_uint64, double_ptr, size_ptr,
    ]
    function.restype = ctypes.c_int
    return function


def _usize(values) -> np.ndarray:
    return np.ascontiguousarray(values, dtype=np.uintp)


def grow(
    monomers,
    conformer_pools,
    *,
    bond_length: float,
    dihedral: float,
    random_rot: bool,
    dist_min: float,
    check_bond_length: bool,
    retry_step: int,
    retry: int,
    rollback: int,
    seed: int | None,
    restart_initial: bool,
):
    """Grow an ordered monomer sequence and return coordinates and diagnostics.

    Monomers and conformer pools are position-aligned. Each pool may have a
    different atom count, which permits copolymer sequences.
    """
    if len(monomers) == 0 or len(monomers) != len(conformer_pools):
        raise ValueError("Rust growth requires one conformer pool per monomer")

    coordinate_parts = []
    coordinate_offsets = [0]
    atom_counts = []
    conformer_counts = []
    heads, tails, head_neighbors, tail_neighbors = [], [], [], []
    atomic_numbers = []
    bond_offsets = [0]
    bond_begin, bond_end, bond_orders = [], [], []

    for index, (monomer, pool) in enumerate(zip(monomers, conformer_pools)):
        mol = monomer.mol
        atom_count = mol.GetNumAtoms()
        if atom_count < 2 or not pool:
            raise ValueError(f"Monomer {index} has no usable coordinates")
        conformers = np.asarray(pool, dtype=np.float64)
        if conformers.ndim != 3 or conformers.shape[1:] != (atom_count, 3):
            raise ValueError(f"Monomer {index} has an invalid conformer-pool shape")
        if not np.isfinite(conformers).all():
            raise ValueError(f"Monomer {index} has non-finite coordinates")

        head = mol.GetIntProp("head_idx")
        tail = mol.GetIntProp("tail_idx")
        head_neighbor = mol.GetIntProp("head_ne_idx")
        tail_neighbor = mol.GetIntProp("tail_ne_idx")
        if not all(0 <= value < atom_count for value in (head, tail, head_neighbor, tail_neighbor)):
            raise ValueError(f"Monomer {index} has invalid linker indices")

        coordinate_parts.append(np.ascontiguousarray(conformers).reshape(-1))
        coordinate_offsets.append(coordinate_offsets[-1] + conformers.size)
        atom_counts.append(atom_count)
        conformer_counts.append(len(conformers))
        heads.append(head)
        tails.append(tail)
        head_neighbors.append(head_neighbor)
        tail_neighbors.append(tail_neighbor)
        atomic_numbers.extend(atom.GetAtomicNum() for atom in mol.GetAtoms())
        for bond in mol.GetBonds():
            bond_begin.append(bond.GetBeginAtomIdx())
            bond_end.append(bond.GetEndAtomIdx())
            bond_orders.append(bond.GetBondTypeAsDouble())
        bond_offsets.append(len(bond_begin))

    coords = np.ascontiguousarray(np.concatenate(coordinate_parts), dtype=np.float64)
    atom_counts_array = _usize(atom_counts)
    conformer_counts_array = _usize(conformer_counts)
    coord_offsets_array = _usize(coordinate_offsets)
    heads_array = _usize(heads)
    tails_array = _usize(tails)
    head_neighbors_array = _usize(head_neighbors)
    tail_neighbors_array = _usize(tail_neighbors)
    atomic_numbers_array = np.ascontiguousarray(atomic_numbers, dtype=np.uint8)
    bond_offsets_array = _usize(bond_offsets)
    bond_begin_array = _usize(bond_begin)
    bond_end_array = _usize(bond_end)
    bond_orders_array = np.ascontiguousarray(bond_orders, dtype=np.float64)
    output = np.empty(sum(atom_counts) * 3, dtype=np.float64)
    stats = np.zeros(3, dtype=np.uintp)

    seed_words = np.random.SeedSequence(seed).generate_state(2, dtype=np.uint32)
    native_seed = int(seed_words[0]) | (int(seed_words[1]) << 32)
    double_ptr = ctypes.POINTER(ctypes.c_double)
    size_ptr = ctypes.POINTER(ctypes.c_size_t)
    u8_ptr = ctypes.POINTER(ctypes.c_uint8)
    status = _core()(
        len(monomers),
        len(coords),
        len(atomic_numbers_array),
        len(bond_begin_array),
        coords.ctypes.data_as(double_ptr),
        coord_offsets_array.ctypes.data_as(size_ptr),
        atom_counts_array.ctypes.data_as(size_ptr),
        conformer_counts_array.ctypes.data_as(size_ptr),
        heads_array.ctypes.data_as(size_ptr),
        tails_array.ctypes.data_as(size_ptr),
        head_neighbors_array.ctypes.data_as(size_ptr),
        tail_neighbors_array.ctypes.data_as(size_ptr),
        atomic_numbers_array.ctypes.data_as(u8_ptr),
        bond_offsets_array.ctypes.data_as(size_ptr),
        bond_begin_array.ctypes.data_as(size_ptr),
        bond_end_array.ctypes.data_as(size_ptr),
        bond_orders_array.ctypes.data_as(double_ptr),
        float(bond_length),
        float(dihedral),
        int(bool(random_rot)),
        float(dist_min),
        int(bool(check_bond_length)),
        int(retry_step),
        int(retry),
        int(rollback),
        int(bool(restart_initial)),
        native_seed,
        output.ctypes.data_as(double_ptr),
        stats.ctypes.data_as(size_ptr),
    )

    attempts, rollbacks, failed_step = (int(value) for value in stats)
    if status == 1:
        return None, attempts, rollbacks, failed_step
    if status != 0:
        raise RuntimeError(f"Rust growth core rejected the input (status {status})")

    coordinates = []
    cursor = 0
    for atom_count in atom_counts:
        end = cursor + atom_count * 3
        coordinates.append(output[cursor:end].reshape(atom_count, 3).copy())
        cursor = end
    return coordinates, attempts, rollbacks, 0
