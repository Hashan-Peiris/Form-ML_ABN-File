import glob
import re
from typing import List, Tuple, Dict

star_line = "**************************************************"
dash_line = "--------------------------------------------------"
eq_line = "=================================================="

# Minimal atomic weights table (u). Extend as needed. Unknown elements default to 0.0
ATOMIC_WEIGHTS: Dict[str, float] = {
    # common elements
    "H": 1.008, "He": 4.002602, "Li": 6.94, "Be": 9.0121831, "B": 10.81, "C": 12.011,
    "N": 14.007, "O": 15.999, "F": 18.998403163, "Ne": 20.1797, "Na": 22.98976928,
    "Mg": 24.305, "Al": 26.9815385, "Si": 28.085, "P": 30.973761998, "S": 32.06,
    "Cl": 35.45, "Ar": 39.948, "K": 39.0983, "Ca": 40.078, "Sc": 44.955908,
    "Ti": 47.867, "V": 50.9415, "Cr": 51.9961, "Mn": 54.938044, "Fe": 55.845,
    "Co": 58.933194, "Ni": 58.6934, "Cu": 63.546, "Zn": 65.38, "Ga": 69.723,
    "Ge": 72.630, "As": 74.921595, "Se": 78.971, "Br": 79.904, "Kr": 83.798,
    "Rb": 85.4678, "Sr": 87.62, "Y": 88.90584, "Zr": 91.224, "Nb": 92.90637,
    "Mo": 95.95, "Tc": 98.0, "Ru": 101.07, "Rh": 102.90550, "Pd": 106.42,
    "Ag": 107.8682, "Cd": 112.414, "In": 114.818, "Sn": 118.710, "Sb": 121.760,
    "Te": 127.60, "I": 126.90447, "Xe": 131.293, "Cs": 132.90545196, "Ba": 137.327,
    "La": 138.90547, "Ce": 140.116, "Pr": 140.90766, "Nd": 144.242, "Sm": 150.36,
    "Eu": 151.964, "Gd": 157.25, "Tb": 158.92535, "Dy": 162.500, "Ho": 164.93033,
    "Er": 167.259, "Tm": 168.93422, "Yb": 173.045, "Lu": 174.9668, "Hf": 178.49,
    "Ta": 180.94788, "W": 183.84, "Re": 186.207, "Os": 190.23, "Ir": 192.217,
    "Pt": 195.084, "Au": 196.966569, "Hg": 200.592, "Tl": 204.38, "Pb": 207.2,
    "Bi": 208.98040
}


def _wrap_three_per_line(values: List[str]) -> str:
    lines: List[str] = []
    for i in range(0, len(values), 3):
        lines.append(" ".join(values[i:i+3]))
    return "\n".join(lines)


def read_poscar(filename: str) -> Tuple[str, float, List[List[float]], List[str], List[int], str, List[List[float]]]:
    with open(filename, "r") as file:
        lines = [ln.rstrip() for ln in file.readlines()]

    if len(lines) < 8:
        raise ValueError(f"POSCAR file {filename} is too short")

    line_idx = 0
    comment = lines[line_idx].strip(); line_idx += 1

    try:
        scale = float(lines[line_idx].strip()); line_idx += 1
    except ValueError:
        # Fallback: treat as 1.0
        scale = 1.0; line_idx += 0

    # Lattice vectors
    cell_vectors = [list(map(float, lines[line_idx + i].split())) for i in range(3)]
    line_idx += 3

    # Atom types list (may be element names or custom like O1, O2)
    atom_types = lines[line_idx].split(); line_idx += 1

    # Atom counts per type
    num_atoms = list(map(int, lines[line_idx].split())); line_idx += 1

    # Handle optional Selective Dynamics line
    coordinate_type = lines[line_idx].strip(); line_idx += 1
    if coordinate_type.lower().startswith("s"):
        # Selective dynamics present; next line is actual coordinate type
        coordinate_type = lines[line_idx].strip(); line_idx += 1

    # Remaining lines are atomic positions (may include T/F flags)
    positional_data_raw: List[List[float]] = []
    for ln in lines[line_idx: line_idx + sum(num_atoms)]:
        parts = ln.split()
        if len(parts) < 3:
            continue
        positional_data_raw.append([float(parts[0]), float(parts[1]), float(parts[2])])

    # Apply scaling to lattice vectors
    scaled_cell_vectors = [[scale * v for v in row] for row in cell_vectors]

    # Convert positions to Cartesian if in fractional/direct
    coord_ch = coordinate_type[:1].lower()
    if coord_ch == "d":
        # fractional -> Cartesian using scaled lattice vectors
        a, b, c = scaled_cell_vectors
        cart_positions: List[List[float]] = []
        for fx, fy, fz in positional_data_raw:
            x = fx * a[0] + fy * b[0] + fz * c[0]
            y = fx * a[1] + fy * b[1] + fz * c[1]
            z = fx * a[2] + fy * b[2] + fz * c[2]
            cart_positions.append([x, y, z])
    else:
        # already Cartesian, but still apply global scale
        cart_positions = [[scale * p for p in row] for row in positional_data_raw]

    return comment, scale, scaled_cell_vectors, atom_types, num_atoms, coordinate_type, cart_positions


def read_outcar_stream(filename: str) -> Tuple[List[List[float]], float, List[float]]:
    forces: List[List[float]] = []
    energy_val: float = None  # type: ignore
    stress_vals: List[float] = []

    # States for forces section
    in_forces = False

    # Regex precompile
    re_forces_hdr = re.compile(r"POSITION\s+TOTAL-FORCE\s*\(eV/Angst\)")
    re_energy = re.compile(r"energy\(sigma->0\)\s*=\s*([\-0-9.Ee+]+)")
    re_energy_alt = re.compile(r"free\s+energy\s+TOTEN\s*=\s*([\-0-9.Ee+]+)")
    re_kb = re.compile(r"in\s+kB(.*)")

    with open(filename, "r") as fh:
        for raw in fh:
            line = raw.rstrip("\n")

            # Energy (keep last occurrence encountered)
            m_e = re_energy.search(line)
            if m_e:
                try:
                    energy_val = float(m_e.group(1))
                except Exception:
                    pass
            else:
                m_e2 = re_energy_alt.search(line)
                if m_e2:
                    try:
                        energy_val = float(m_e2.group(1))
                    except Exception:
                        pass

            # Stress line (last one wins)
            m_kb = re_kb.search(line)
            if m_kb:
                tail = m_kb.group(1).strip()
                try:
                    cand = list(map(float, tail.split()))
                    if len(cand) >= 6:
                        stress_vals = cand[:6]
                except Exception:
                    pass

            # Forces block detection and extraction
            if not in_forces:
                if re_forces_hdr.search(line):
                    # The next line is a separator/header; actual data starts after that
                    in_forces = True
                    # Skip the very next line (header underline) by reading one more
                    next(fh, None)
                    continue
            else:
                # End of forces section marked by a line containing 'total drift:'
                if "total drift:" in line:
                    in_forces = False
                    continue
                if not line.strip():
                    continue
                parts = line.split()
                if len(parts) >= 6:
                    # columns: x y z fx fy fz
                    try:
                        fx, fy, fz = float(parts[3]), float(parts[4]), float(parts[5])
                        forces.append([fx, fy, fz])
                    except Exception:
                        pass

    if energy_val is None:
        raise RuntimeError(f"Failed to parse total energy from {filename}")
    if len(stress_vals) != 6:
        raise RuntimeError(f"Failed to parse 6-component stress from {filename}")

    return forces, energy_val, stress_vals


def write_header(outf, num_configs: int, all_types_in_order: List[str], max_atoms_per_system: int,
                 max_atoms_per_atom_type_scalar: int) -> None:
    # 1.0 Version
    outf.write("1.0 Version\n")
    outf.write(f"{star_line}\n")

    # The number of configurations
    outf.write("The number of configurations\n")
    outf.write(f"{dash_line}\n")
    outf.write(f"{num_configs}\n")
    outf.write(f"{star_line}\n")

    # The maximum number of atom type (total unique types across all structures)
    outf.write("The maximum number of atom type\n")
    outf.write(f"{dash_line}\n")
    outf.write(f"{len(all_types_in_order)}\n")
    outf.write(f"{star_line}\n")

    # The atom types in the data file (max 3 per line)
    outf.write("The atom types in the data file\n")
    outf.write(f"{dash_line}\n")
    outf.write(f"{_wrap_three_per_line(all_types_in_order)}\n")
    outf.write(f"{star_line}\n")

    # The maximum number of atoms per system
    outf.write("The maximum number of atoms per system\n")
    outf.write(f"{dash_line}\n")
    outf.write(f"{max_atoms_per_system}\n")
    outf.write(f"{star_line}\n")

    # The maximum number of atoms per atom type (single scalar)
    outf.write("The maximum number of atoms per atom type\n")
    outf.write(f"{dash_line}\n")
    outf.write(f"{max_atoms_per_atom_type_scalar}\n")
    outf.write(f"{star_line}\n")

    # Reference atomic energy (eV) — element dependent, here zeros as placeholder
    outf.write("Reference atomic energy (eV)\n")
    outf.write(f"{dash_line}\n")
    ref_energies = [f"{0.0}" for _ in all_types_in_order]
    outf.write(f"{_wrap_three_per_line(ref_energies)}\n")
    outf.write(f"{star_line}\n")

    # Atomic mass — element dependent, in u
    outf.write("Atomic mass\n")
    outf.write(f"{dash_line}\n")
    masses = [f" {ATOMIC_WEIGHTS.get(t, 0.0)}" for t in all_types_in_order]
    outf.write(f"{_wrap_three_per_line(masses)}\n")
    outf.write(f"{star_line}\n")

    # The numbers of basis sets per atom type — element dependent; for select mode dummy 1s
    outf.write("The numbers of basis sets per atom type \n")
    outf.write(f"{dash_line}\n")
    basis_counts = ["1" for _ in all_types_in_order]
    outf.write(f"{_wrap_three_per_line(basis_counts)}\n")
    outf.write(f"{star_line}\n")

    # Basis set for X — dummy block per species with one entry "1 1"
    for t in all_types_in_order:
        outf.write(f"Basis set for {t} \n")
        outf.write(f"{dash_line}\n")
        outf.write("1 1\n")
        outf.write(f"{star_line}\n")


def write_configuration(outf, index_one_based: int, system_name: str, cell_vectors: List[List[float]],
                        atom_types: List[str], num_atoms: List[int], positions_cart: List[List[float]],
                        total_energy: float, forces: List[List[float]], stress_kbar: List[float]) -> None:
    outf.write(f"Configuration num. {index_one_based}\n")
    outf.write(f"{eq_line}\n")

    # System name from POSCAR comment
    outf.write("System name\n")
    outf.write(f"{dash_line}\n")
    outf.write(f"{system_name[:40]}\n")
    outf.write(f"{eq_line}\n")

    # Number of atom types in this structure
    outf.write("The number of atom types\n")
    outf.write(f"{dash_line}\n")
    outf.write(f"{len(atom_types)}\n")
    outf.write(f"{eq_line}\n")

    # Number of atoms total
    total_atoms = sum(num_atoms)
    outf.write("The number of atoms\n")
    outf.write(f"{dash_line}\n")
    outf.write(f"{total_atoms}\n")
    outf.write(f"{star_line}\n")

    # Atom types and atom numbers (per line)
    outf.write("Atom types and atom numbers\n")
    outf.write(f"{dash_line}\n")
    for t, n in zip(atom_types, num_atoms):
        outf.write(f" {t} {n}\n")
    outf.write(f"{eq_line}\n")

    # Primitive lattice vectors (ang.)
    outf.write("Primitive lattice vectors (ang.)\n")
    outf.write(f"{dash_line}\n")
    for row in cell_vectors:
        outf.write(" ".join(str(v) for v in row) + "\n")
    outf.write(f"{eq_line}\n")

    # Atomic positions (ang.) — Cartesian
    outf.write("Atomic positions (ang.)\n")
    outf.write(f"{dash_line}\n")
    for pos in positions_cart:
        outf.write(" ".join(str(v) for v in pos) + "\n")
    outf.write(f"{eq_line}\n")

    # Total energy (eV)
    outf.write("Total energy (eV)\n")
    outf.write(f"{dash_line}\n")
    outf.write(f"{total_energy}\n")
    outf.write(f"{eq_line}\n")

    # Forces (eV ang.^-1)
    outf.write("Forces (eV ang.^-1)\n")
    outf.write(f"{dash_line}\n")
    for fx, fy, fz in forces:
        outf.write(f"  {fx} {fy} {fz}\n")
    outf.write(f"{eq_line}\n")

    # Stress (kbar)
    outf.write("Stress (kbar)\n")
    outf.write(f"{dash_line}\n")
    outf.write("XX YY ZZ \n")
    outf.write(f"{dash_line} \n")
    outf.write(f" {' '.join(str(x) for x in stress_kbar[:3])} \n")
    outf.write(f"{dash_line}\n")
    outf.write("XY YZ ZX \n")
    outf.write(f"{dash_line} \n")
    outf.write(f"{' '.join(str(x) for x in stress_kbar[3:6])} \n")
    outf.write(f"{star_line}\n")


def main() -> None:
    poscar_files = sorted(glob.glob("POSCAR_*"))
    outcar_files = sorted(glob.glob("OUTCAR_*"))
    assert len(poscar_files) == len(outcar_files), "Mismatch in number of POSCAR and OUTCAR files."
    assert len(poscar_files) > 0, "No POSCAR_*/OUTCAR_* pairs found."

    # First pass: collect header info from POSCAR only (fast)
    all_types_in_order: List[str] = []
    seen_types: set = set()
    max_atoms_per_system = 0
    max_atoms_per_atom_type_scalar = 0

    # For per-configuration typing, we need to honor the order in POSCAR (as given in line 6)
    for poscar_file in poscar_files:
        comment, scale, cell_vectors, elements_data, num_atoms, coordinate_type, positions = read_poscar(poscar_file)
        # Track union of types preserving order of first appearance
        for t in elements_data:
            if t not in seen_types:
                seen_types.add(t)
                all_types_in_order.append(t)
        # Max atoms per system
        total_atoms = sum(num_atoms)
        if total_atoms > max_atoms_per_system:
            max_atoms_per_system = total_atoms
        # Max atoms per atom type (single scalar across all types and structures)
        local_max = max(num_atoms) if num_atoms else 0
        if local_max > max_atoms_per_atom_type_scalar:
            max_atoms_per_atom_type_scalar = local_max

    # Write header and configurations streaming to disk
    with open("ML_ABN", "w") as outf:
        write_header(outf, num_configs=len(poscar_files), all_types_in_order=all_types_in_order,
                     max_atoms_per_system=max_atoms_per_system,
                     max_atoms_per_atom_type_scalar=max_atoms_per_atom_type_scalar)

        # Second pass: write configurations; parse OUTCAR on the fly
        for idx, (poscar_file, outcar_file) in enumerate(zip(poscar_files, outcar_files), start=1):
            comment, scale, cell_vectors, elements_data, num_atoms, coordinate_type, positions_cart = read_poscar(poscar_file)
            forces, energy, stress = read_outcar_stream(outcar_file)
            write_configuration(outf, idx, comment, cell_vectors, elements_data, num_atoms,
                                positions_cart, energy, forces, stress)


if __name__ == "__main__":
    main()