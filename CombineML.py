import glob
import re

star_line = "**************************************************"
dash_line = "--------------------------------------------------"
eq_line   = "=================================================="

def read_poscar(filename):
    with open(filename, "r") as file:
        lines = file.readlines()
    
    comment = lines[0].strip()
    scale = float(lines[1].strip())
    cell_vectors = [list(map(float, line.strip().split())) for line in lines[2:5]]
    atom_types = lines[5].strip().split()
    num_atoms = list(map(int, lines[6].strip().split()))
    coordinate_type = lines[7].strip()
    atomic_positions = [list(map(float, line.strip().split())) for line in lines[8:]]
    
    return comment, scale, cell_vectors, atom_types, num_atoms, coordinate_type, atomic_positions

def read_outcar(filename):
    with open(filename, "r") as file:
        content = file.read()
    
    forces_section_start = content.index(re.search(r"POSITION\s+TOTAL-FORCE \(eV/Angst\)", content).group(0))
    content_lines = content[forces_section_start:].split("\n")
    end_line_keyword = "total drift:"
    forces_section_end = next(i for i, line in enumerate(content_lines[2:]) if end_line_keyword in line.strip())
    
    forces_section = content_lines[2:forces_section_end+1]
    forces_data = []
    for line in forces_section:
        if not line.strip():
            continue
        line_data = list(map(float, line.split()[3:6])) # Get data from positions 3 to 5 (Force X, Y, Z)
        forces_data.append(line_data)
    
    energy_line = re.search(r"energy\(sigma->0\)\s+=\s+(.*?)\n", content).group(1)
    energy_data = float(energy_line.strip())
    
    stress_line = re.search(r"in kB(.*?)\n", content).group(1)
    stress_data = list(map(float, stress_line.strip().split()))
    
    return forces_data, energy_data, stress_data

def form_ml_abn_file(header_data, config_data):
    version_text = f"1.0 Version\n{star_line}"
    num_configs_text = f"The number of configurations\n{dash_line}\n{header_data['num_configs']}\n{star_line}"
    max_num_atomtypes_text = f"The maximum number of atom type\n{dash_line}\n{header_data['max_num_atomtypes']}\n{star_line}"
    atomtypes_text = f"The atom types in the data file\n{dash_line}\n{' '.join(header_data['atomtypes'])}\n{star_line}"
    max_num_atoms_text = f"The maximum number of atoms per system\n{dash_line}\n{header_data['max_num_atoms']}\n{star_line}"
    max_atoms_per_atomtype_text = f"The maximum number of atoms per atom type\n{dash_line}\n" + "\n".join(f"{max_cnt}" for max_cnt in header_data['max_atoms_per_atomtype'].values()) + f"\n{star_line}"
    
    config_texts = []
    for i, config in enumerate(config_data, start=1):
        config_texts.append(f"Configuration num. {i}")
        config_texts.append(f"{eq_line}")
        config_texts.append(f'System name\n{dash_line}\nSystem\n{eq_line}')
        config_texts.append(f"The number of atom types\n{dash_line}\n{len(set(config['elements_data']))}\n{eq_line}")
        config_texts.append(f"The number of atoms\n{dash_line}\n{sum(config['num_atoms'])}\n{star_line}")
        config_texts.append(f'Atom types and atom numbers\n{dash_line}')
        
        element_counts = {}
        for element, count in zip(config['elements_data'], config['num_atoms']):
            if element in element_counts:
                element_counts[element] += count
            else:
                element_counts[element] = count

        for element, count in element_counts.items():
            config_texts.append(f" {element} {count}")
        
        config_texts.append(f"{eq_line}\nPrimitive lattice vectors (ang.)\n{dash_line}\n" + "\n".join(" ".join(map(str, line)) for line in config['cell_vectors']) + f"\n{eq_line}")
        config_texts.append(f"Atomic positions (ang.)\n{dash_line}\n" + "\n".join(" ".join(map(str, line)) for line in config['positional_data']) + f"\n{eq_line}")
        config_texts.append(f"Total energy (eV)\n{dash_line}\n{config['energy_data']}\n{eq_line}")
        config_texts.append(f"Forces (eV ang.^-1)\n{dash_line}")
        for force in config['force_data']:
            config_texts.append(f"  {' '.join(map(str, force))}")
        config_texts.append(f"{eq_line}\nStress (kbar)\n{dash_line}\nXX YY ZZ \n{dash_line} \n {' '.join(map(str, config['stress_data'][:3]))} \n{dash_line}")
        config_texts.append(f"XY YZ ZX \n{dash_line} \n{' '.join(map(str, config['stress_data'][3:]))} \n{star_line}")
    
    ref_atomic_energy_text = f"Reference atomic energy (eV)\n{dash_line}\n0.0\n{star_line}"
    atomic_masses = {'Au': 197.0}  # Placeholder dictionary, replace with actual atomic masses
    atomic_mass_text = f"Atomic mass\n{dash_line}\n" + "\n".join(f" {mass}" for mass in atomic_masses.values()) + f"\n{star_line}"
    basis_sets = f"The numbers of basis sets per atom type \n{dash_line}\n1\n{star_line}"
    basis_set_num = f"Basis set for Au \n{dash_line}\n1 1\n{star_line}"

    ml_abn_content = "\n".join([version_text, num_configs_text, max_num_atomtypes_text, atomtypes_text, max_num_atoms_text,
                                max_atoms_per_atomtype_text, ref_atomic_energy_text, atomic_mass_text, basis_sets, basis_set_num] + config_texts)
    
    return ml_abn_content

def main():
    poscar_files = glob.glob("POSCAR_*")
    outcar_files = glob.glob("OUTCAR_*")
    poscar_files.sort()
    outcar_files.sort()
    assert len(poscar_files) == len(outcar_files), "Mismatch in number of POSCAR and OUTCAR files."

    header_data = {}
    config_data = []
    
    max_num_atoms = -float("inf")
    max_num_atomtypes = -float("inf")
    atomtypes = set()
    max_atoms_per_type = {}

    for poscar_file, outcar_file in zip(poscar_files, outcar_files):
        comment, scale, cell_vectors, elements_data, num_atoms, coordinate_type, positional_data = read_poscar(poscar_file)
        force_data, energy_data, stress_data = read_outcar(outcar_file)
        
        max_num_atoms = max(max_num_atoms, sum(num_atoms))
        max_num_atomtypes = max(max_num_atomtypes, len(set(elements_data)))
        atomtypes.update(elements_data)
        
        for element, count in zip(elements_data, num_atoms):
            if element not in max_atoms_per_type:
                max_atoms_per_type[element] = count
            else:
                max_atoms_per_type[element] = max(max_atoms_per_type[element], count)
        
        config_data.append({
            'comment': comment,
            'scale': scale,
            'cell_vectors': cell_vectors,
            'elements_data': elements_data,
            'num_atoms': num_atoms,
            'coordinate_type': coordinate_type,
            'positional_data': positional_data,
            'force_data': force_data,
            'energy_data': energy_data,
            'stress_data': stress_data
        })

    header_data['num_configs'] = len(poscar_files)
    header_data['max_num_atoms'] = max_num_atoms
    header_data['max_num_atomtypes'] = max_num_atomtypes
    header_data['atomtypes'] = list(atomtypes)
    header_data['max_atoms_per_atomtype'] = max_atoms_per_type

    ml_abn_content = form_ml_abn_file(header_data, config_data)

    with open("ML_ABN", "w") as file:
        file.write(ml_abn_content)

if __name__ == "__main__":
    main()