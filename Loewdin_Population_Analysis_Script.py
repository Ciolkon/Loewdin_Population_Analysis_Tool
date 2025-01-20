import pandas as pd
import matplotlib.pyplot as plt
import os

# Change working directory to the location of the ORCA output file
os.chdir(r"C:\Users\nickc\Documents\Data\Compiled_data\FeDimer\Calculations\KLFeOPeroxo_SCF_UNO_SVP")
print(os.getcwd())  # Prints the current working directory

output_file_path = r"C:\Users\nickc\Documents\Data\Compiled_data\FeDimer\Calculations\KLFeOPeroxo_SCF_UNO_SVP\KLFeOPeroxo_SCF_UNO.out"

if os.path.exists(output_file_path):
    print("File found!")
else:
    print("File not found! Check the path.")

def parse_loewdin_population(output_lines, mo_filter):
    """
    Parses the Loewdin population section of the ORCA output file, handling pagination.
    Retains lines of dashes as page breaks after relevant MOs are recorded.
    """
    loewdin_data = []
    capture = False
    mo_headers = []
    header_indices = {}
    relevant_mo_detected = False  # Flag to track when relevant MOs are encountered

    for line in output_lines:
        stripped_line = line.strip()

        # Detect start of the Loewdin population analysis section
        if "LOEWDIN REDUCED ORBITAL POPULATIONS PER UNO" in line.upper():
            capture = True
            continue

        # Retain lines of dashes as page breaks, but only after relevant MOs are detected
        if capture and all(char == "-" for char in stripped_line):
            if relevant_mo_detected:  # Add page break only if relevant MOs have been processed
                loewdin_data.append(("PAGE_BREAK", stripped_line, []))
            continue

        # Detect and parse MO headers (rows that contain only MO numbers)
        if capture and all(x.isdigit() for x in stripped_line.split()):
            mo_headers = [int(x) for x in stripped_line.split()]
            header_indices = {mo: idx for idx, mo in enumerate(mo_headers) if mo in mo_filter}
            if header_indices:  # If there are MOs in `mo_filter`, set the relevant flag
                relevant_mo_detected = True
            continue

        # Parse data rows if MO headers have been found
        if capture and header_indices:
            parts = stripped_line.split()
            try:
                # Extract the populations for filtered MOs
                populations = [float(parts[3 + idx]) for mo, idx in header_indices.items()]
                # Extract atom and orbital information
                atom_number, atom_type, orbital = parts[:3]
                loewdin_data.append((f"{atom_number} {atom_type}", orbital, populations))
            except (ValueError, IndexError):
                # Skip lines that cannot be parsed
                continue

        # Stop capture after processing all specified MOs
        if capture and not header_indices and mo_headers and max(mo_headers) >= max(mo_filter):
            break

    return loewdin_data

def write_loewdin_parsed_to_file(loewdin_parsed, output_filename):
    """
    Writes the parsed Loewdin data to a text file.

    Parameters:
    - loewdin_parsed (list of tuples): Parsed Loewdin population data.
    - output_filename (str): Name of the output text file.
    """
    with open(output_filename, 'w') as file:
        file.write("Atom Info, Orbital, Populations\n")
        for atom_info, orbital, populations in loewdin_parsed:
            populations_str = ", ".join(map(str, populations))
            file.write(f"{atom_info}, {orbital}, {populations_str}\n")
    print(f"Loewdin parsed data written to {output_filename}")

def process_orca_output(file_path, atoms_to_analyze, mo_filter):
    """
    Processes the ORCA output file to sum the provided percentages
    for specified atoms and orbitals across MOs.
    """
    # Read the content of the file
    with open(file_path, 'r') as file:
        output_content = file.readlines()

    # Parse the Loewdin population data
    loewdin_parsed = parse_loewdin_population(output_content, mo_filter)
    write_loewdin_parsed_to_file(loewdin_parsed, "parsing_debug.txt")
    
    # Initialize a dictionary to store summed contributions
    summary = {atom: {mo: 0 for mo in mo_filter} for atom in atoms_to_analyze.keys()}
    individual_contributions = {}  # For individual atom contributions

    # Initialize a row offset to track MO alignment across rows
    row_offset = 0
    last_mo_index = 0  # Tracks the index of the last MO processed

    # Loop through the parsed data and aggregate contributions
    for entry in loewdin_parsed:
        # Detect and handle page breaks
        if entry[0] == "PAGE_BREAK":
            row_offset = last_mo_index  # Adjust the row offset for the new page
            last_mo_index = 0  # Reset the last_mo_index for the next page
            continue

        # Validate the data structure and skip invalid entries
        if len(entry) < 3 or not isinstance(entry[2], list):
            print(f"Skipping invalid entry: {entry}")
            continue
          
        # Process the normal atomic data rows
        try:
            atom_info, orbital, populations = entry
            atom_number, atom_type = atom_info.split()
            atom_number = int(atom_number)  # Ensure atom_number is an integer
        except ValueError as e:
            print(f"Error parsing entry {entry}: {e}")
            continue  # Skip invalid rows

        # Debug: Log the current row offset
        print(f"Atom: {atom_info}, Orbital: {orbital}, Populations: {populations}, Row Offset: {row_offset}")

        # Check if the atom type is in the filters
        if atom_type in atoms_to_analyze:
            allowed_orbitals = atoms_to_analyze[atom_type]["orbitals"]
            specific_atoms = atoms_to_analyze[atom_type]["atoms"]
            
            # Check if the orbital matches the filters
            if any(orbital_filter in orbital for orbital_filter in allowed_orbitals):
                # Check if specific atom numbers are provided, or include all atoms of that type
                if not specific_atoms or atom_number in specific_atoms:
                    # Iterate through populations and match them to mo_filter using row_offset
                    for i, pop in enumerate(populations):
                        mo_index = (row_offset + i) % len(mo_filter)  # Adjust MO index with row_offset
                        mo = mo_filter[mo_index]
                        
                        # Add contributions to the summary
                        summary[atom_type][mo] += pop
                        
                        # Track individual contributions if specific atoms are specified
                        if specific_atoms:
                            individual_contributions.setdefault(atom_number, {}).setdefault(mo, 0)
                            individual_contributions[atom_number][mo] += pop
                        
                        last_mo_index = mo_index + 1  # Update last_mo_index

    # Convert summary to DataFrame for visualization
    summary_df = pd.DataFrame(summary).T
    summary_df.columns = [f"MO{i}" for i in mo_filter]

    # Convert individual contributions to DataFrame if applicable
    if individual_contributions:
        individual_df = pd.DataFrame.from_dict(individual_contributions, orient='index')
        individual_df.columns = [f"MO{i}" for i in mo_filter]
    else:
        individual_df = None

    print(summary_df)
    print(individual_df)
    return summary_df, individual_df

# Inputs
atoms_to_analyze = {
    "Fe": {"orbitals": ["d"], "atoms": []},   # Element, Atomic orbitals, atom index number
    "O": {"orbitals": ["s", "p"], "atoms": []},  
    "N": {"orbitals": ["s", "p"], "atoms": []},  
    "P": {"orbitals": ["s", "p"], "atoms": []}   
}

mos_to_analyze = [520, 521, 522, 523, 524, 525, 526, 527, 528, 529]  # Specify MO indices

# Process the file
summary_df, individual_df = process_orca_output(output_file_path, atoms_to_analyze, mos_to_analyze)

# Plot the overall results as a line graph
plt.figure(figsize=(10, 6))
for atom_type in summary_df.index:
    plt.plot(summary_df.columns, summary_df.loc[atom_type], marker='o', label=atom_type)

# If individual contributions exist, plot them as separate lines
if individual_df is not None:
    for atom_number in individual_df.index:
        # Find the corresponding atom_type for the atom_number
        matched_atom_type = None
        for atom_type, atom_info in atoms_to_analyze.items():
            if atom_number in atom_info["atoms"]:
                matched_atom_type = atom_type
                break

        if matched_atom_type:
            plt.plot(
                individual_df.columns,
                individual_df.loc[atom_number],
                '--',
                label=f"Atom {atom_number} ({matched_atom_type})"
            )
        else:
            print(f"Warning: Atom {atom_number} not found in atoms_to_analyze mapping.")

plt.title('Percentage Contributions to Molecular Orbitals (MOs)')
plt.ylabel('Percentage Contribution (%)')
plt.xlabel('MOs')
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.legend(title='Atoms/Atom Types', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.show()
