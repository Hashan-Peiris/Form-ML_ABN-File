import glob
import re

def read_poscar(filename):
    # Reading the POSCAR file
    with open(filename, "r") as file:
        lines = file.readlines()
    
    # The first line is the comment
    comment = lines[0].strip()
    
    # The second line is the scale
    scale = float(lines[1].strip())
    
    # The next three lines are the cell vectors
    cell_vectors = [list(map(float, line.strip().split())) for line in lines[2:5]]
    
    # The next line is the atom types
    atom_types = lines[5].strip().split()

    # The next line is the number of atoms
    num_atoms = list(map(int, lines[6].strip().split()))

    # The next line is the coordinate type
    coordinate_type = lines[7].strip()

    # Remaining lines are the atomic position data
    atomic_positions = [list(map(float, line.strip().split())) for line in lines[8:]]

    return comment, scale, cell_vectors, atom_types, num_atoms, coordinate_type, atomic_positions


def read_outcar(filename):
    # Reading the OUTCAR file
    with open(filename, "r") as file:
        content = file.read()

    # Finding and extracting forces data
    forces_section_start = content.index(re.search(r"POSITION\s+TOTAL-FORCE \(eV/Angst\)", content).group(0))
    content_lines = content[forces_section_start:].split("\n")
    end_line_keyword = "total drift:"
    forces_section_end = next(i for i, line in enumerate(content_lines[2:]) if end_line_keyword in line.strip())

    forces_section = content_lines[2:forces_section_end+1]  # Start from 2 to skip the two lines
    forces_data = []
    for line in forces_section:
        if not line.strip():  # Added this line
            continue          # And this one
        line_data = list(map(float, line.split()[3:6]))  # Get data from positions 3 to 5 (Force X, Y, Z)
        forces_data.append(line_data)
    
    print(f"Forces in file {filename}:")
    print(forces_data[:2])  # First two lines of forces
    print(forces_data[-2:])  # Last two lines of forces
    print("------")
    
    # Finding and extracting energy data
    energy_line = re.search(r"energy\(sigma->0\) =(.*?)\n", content).group(1)
    energy_data = float(energy_line.strip())
    
    # Finding and extracting stress data
    stress_line = re.search(r"in kB(.*?)\n", content).group(1)
    stress_data = list(map(float, stress_line.strip().split()))
    
    return forces_data, energy_data, stress_data

# Placeholder function to form ML_ABN file lines
def form_ml_abn_file(header_data, config_data):
    # Formatting the header
    version_text = " 1.0 Version\n"
    num_configs_text = f" The number of configurations\n{header_data['num_configs']}\n"
    max_num_atoms_text = f" The maximum number of atoms per system\n{header_data['max_num_atoms']}\n"
    max_num_atomtypes_text = f" The maximum number of atom type\n{header_data['max_num_atomtypes']}\n"
    atomtypes_text = f" The atom types in the data file\n{' '.join(header_data['atomtypes'])}\n"
    
    # Formatting configs
    config_texts = []
    for i, config in enumerate(config_data, start=1):   # Numbering starts from 1
        config_texts.append(f" Configuration num. {i}")
        config_texts.append(' System name\nSystem')  # Not sure what the System Name should be
        config_texts.append(f" The number of atom types\n{len(set(config['elements_data']))}")  # Assuming elements data contains actual element names
        config_texts.append(f" The number of atoms\n{len(config['elements_data'])}")  # Assuming elements data contains actual element names
        config_texts.append(' Atom types and atom numbers')
        for element in set(config['elements_data']):
            config_texts.append(f" {element} {config['elements_data'].count(element)}")  # Number of each type of atom
        #CTIFOR value will be ignored and not written
        #config_texts.append(" CTIFOR (optional)\n100")  
        config_texts.append(" Primitive lattice vectors (ang.)\n" + "\n".join(" ".join(map(str, line)) for line in config['cell_vectors']))
        config_texts.append(" Atomic positions (ang.)\n" + "\n".join(" ".join(map(str, line)) for line in config['positional_data']))
        config_texts.append(f" Total energy (eV)\n{config['energy_data']}")
        config_texts.append(" Forces (eV ang.^-1)")
        for force in config['force_data']:
            config_texts.append(f"  {' '.join(map(str, force))}")
        config_texts.append(f" Stress (kbar)\n{' '.join(map(str, config['stress_data'][:3]))}")  # Joining first three stress data values
        config_texts.append(f"{' '.join(map(str, config['stress_data'][3:]))}")  # Joining last three stress data values


    # Joining the header and formatted configs
    ml_abn_content = "\n".join([version_text, num_configs_text, max_num_atoms_text, max_num_atomtypes_text, atomtypes_text] + config_texts)
    
    return ml_abn_content 

def main():
    # Loading all the POSCAR and OUTCAR files
    poscar_files = glob.glob("POSCAR_*")
    outcar_files = glob.glob("OUTCAR_*")

    # Sort the files to keep the configurations in order
    poscar_files.sort()
    outcar_files.sort()

    # Validate if POSCAR and OUTCAR files are of equal number
    assert len(poscar_files) == len(outcar_files), "Mismatch in number of POSCAR and OUTCAR files."

    header_data = {}
    config_data = []
    
    # Temporary variables to store max number of atoms and atom types
    max_num_atoms = -float("inf")
    max_num_atomtypes = -float("inf")
    atomtypes = set()

# Process each POSCAR and OUTCAR file pair
    for poscar_file, outcar_file in zip(poscar_files, outcar_files):
        comment, scale, cell_vectors, elements_data, num_atoms, coordinate_type, positional_data = read_poscar(poscar_file)
        force_data, energy_data, stress_data = read_outcar(outcar_file)
        
        # Update the max_num_atoms and max_num_atomtypes if needed
        max_num_atoms = max(max_num_atoms, len(elements_data))
        num_atomtypes = len(set(elements_data))
        max_num_atomtypes = max(max_num_atomtypes, num_atomtypes)
        
        # Add new atom types to the set
        atomtypes.update(elements_data)

        # Gathering the data into a single configuration data
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

        # Setting the header data
        header_data['num_configs'] = len(poscar_files)  # Number of configurations is equal to number of POSCAR files
        header_data['max_num_atoms'] = max_num_atoms
        header_data['max_num_atomtypes'] = max_num_atomtypes
        header_data['atomtypes'] = list(atomtypes)

    # Form the ML_ABN file content
    ml_abn_content = form_ml_abn_file(header_data, config_data)

    # Writing to ML_ABN file
    with open("ML_ABN", "w") as file:
        file.write(ml_abn_content)

if __name__ == "__main__":
    main()
