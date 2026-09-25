"""

Main.py

================================================================================
THESIS TITLE : Analisis in silico del efecto de mutaciones en la interaccion de Thanatin con LptA 
AUTHOR       : Alexandro Mauricio Rizzo Gomez
MASTER       : Master in Bioinformatics
INSTITUTION  : Universidad Internacional de Valencia
SUPERVISORS  : Miquel Sendra
CREATED      : 2026-07-03
LAST UPDATE  : 2026-09-24
================================================================================
FILE PURPOSE:
    This code is part of a Master Thesis for a course in Bioinformatics at the International University of Valencia. 
    This code is written to support the generation of multiple yaml files for the analysis with OSPREY in an automated form. 

================================================================================

"""
#Imports the required libraries for the correct functioning of the program
import csv
import math
import os
import re
from collections import defaultdict

#An amino acid translation dictionary for later use in the code. 
amino_acid_map = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
    'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
    'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
    'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
}

#Some basic variables are defined for the code to work. 
CSV_FILE = "1snp_mutants.csv"
CHAIN_ID = "A"
PDB_OFFSET = 0


def parse_pdb_atoms(yaml_text):
    #This code extracts the information of the coordinates of the atoms of the pdb file into a dictionary fur further processing in the code. 
    #Through the defaultdic function, and lambda the sections are being filled automatically. 
    pdb_data = defaultdict(lambda: defaultdict(lambda: {'res_name': '', 'coords': []}))

    #The code process the yaml file line by line, to extract the data to the dataframe.
    for line in yaml_text.splitlines():
        line_str = line.strip()
        
        # the if statement selects lines with coordinates for atoms
        if line_str.startswith("ATOM") and len(line_str) >= 54:
            
            #the code below extract the data from the lines according to the standard pdb file format. 
            res_name = line_str[17:20].strip()
            chain = line_str[21].strip()
            res_num_str = line_str[22:26].strip()

            #Code for extraction of coordinates data. 
            if res_num_str.lstrip("-").isdigit():
                res_num = int(res_num_str)
                x = float(line_str[30:38].strip())
                y = float(line_str[38:46].strip())
                z = float(line_str[46:54].strip())

                #Data extracte is stored in the pdb_data 
                pdb_data[chain][res_num]['res_name'] = res_name
                pdb_data[chain][res_num]['coords'].append((x, y, z))

    return pdb_data


def find_3d_pocket_neighbors(target_position, target_chain, pdb_data, cutoff_radius):
    #The objective of this piece of coding is to select residues in the proximity of specific residue. 
    
    #It gets the coordinates of a specific residue from the pdb_data dataframte and gets them, for the residue selected (target_position). 
    #it has a safe option to avoid crash, returning nothing if there was a missing residue. 
    target_coords = pdb_data[target_chain].get(target_position, {}).get('coords', [])
    if not target_coords:
        return []

    #A empty list of neihbors residues is created to be filled. 
    neighbors = []
    
    #With a loop, it scans the pdb data for the neigbours, with a bypass system to prevent doing calculations on itself. 
    for chain, residues in pdb_data.items():
        for res_num, res_info in residues.items():
            if chain == target_chain and res_num == target_position:
                continue

            #THe code wth the coordinates of a atom from the target residue (t-xyz) the distance to a neighbour atom candidate (n-xyz) to obtain the distance.
            #If the atom falls in the cutoff radious that the user has selected, the residue is selected as neighbor. And the loop is stop through the break. 
            is_neighbor = False
            for tx, ty, tz in target_coords:
                for nx, ny, nz in res_info['coords']:
                    dist = math.sqrt((tx - nx) ** 2 + (ty - ny) ** 2 + (tz - nz) ** 2)
                    if dist <= cutoff_radius:
                        is_neighbor = True
                        break
                if is_neighbor:
                    break

            #Finally data is appended into the list created at the start of the function. 
            if is_neighbor:
                neighbors.append((chain, res_num, res_info['res_name']))

    return sorted(neighbors, key=lambda x: (x[0], x[1]))


def format_residue_config(chain, res_num, aa_type, mut_list=None):
    #This code builds a formal bloc for the yaml file that will eventually be generated. 
    
    #If the mutation is a target mtuation selected earlier on the workflow, the code sorts alphabetically to create the YAML array of mutations. 
    #If it is a neighbour it assigns a bracket into the YAML file. 
    mut_str = f"[{', '.join(sorted(mut_list))}]" if mut_list else "[]"
    return (
        f"    - mutability: {mut_str}\n"
        f"      flexibility:\n"
        f"        is_flexible: true\n"
        f"        include_structure_rotamer: true\n"
        f"        use_continuous: true\n"
        f"      identity:\n"
        f"        chain: {chain}\n"
        f"        res_num: {res_num}\n"
        f"        aa_type: {aa_type}\n"
    )


def build_3d_flexible_1snp_yaml_files(base_yaml_path, cutoff, csv_path=CSV_FILE, chain_id=CHAIN_ID, offset=PDB_OFFSET):
    
    
    #Error handles for not found files. 
    if not os.path.exists(base_yaml_path):
        print(f"Error: '{base_yaml_path}' not found.")
        return
    if not os.path.exists(csv_path):
        print(f"Error: '{csv_path}' not found. Please execute Main.py first.")
        return

    #This code extracts the name of the file for further processing later on
    yaml_prefix = os.path.splitext(os.path.basename(base_yaml_path))[0]

    with open(base_yaml_path, "r") as f:
        raw_yaml = f.read()

    pdb_data = parse_pdb_atoms(raw_yaml)
    pdb_wt_map = {res: info['res_name'] for res, info in pdb_data[chain_id].items()}
    

    #Cleans any posible old residues: or scan: blocks in the yaml file.
    clean_yaml = re.sub(
        r'scan:\s*\n(?:\s*residues:\s*\n)?(?:\s*scan:\s*\n)?\s*residues:\s*\n(?:\s*-\s*identity:\s*\n\s*chain:\s*\S+\s*\n\s*res_num:\s*\d+\s*\n\s*aa_type:\s*\S+\s*\n)+',
        '',
        raw_yaml
    )
    clean_yaml = re.sub(r'scan:.*?(?=\nprotein:|\nligand:|\Z)', '', clean_yaml, flags=re.DOTALL).strip()

    # This code loads the mutations from the 1snp_mutants.csv. There is a fallback option to provide a safe net instead of provoking a crash
    targets = defaultdict(lambda: {'muts': set(), 'fallback_wt': 'SER'})
    with open(csv_path, mode="r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            raw_pos = str(row.get("position", "")).strip()
            mut_1L = str(row.get("mutant", "")).strip().upper()
            orig_1L = str(row.get("original", "")).strip().upper()

            #Coding for dealing with empty rows or stop codons. 
            if not raw_pos or not mut_1L or "*" in mut_1L:
                continue

            #Coding to deal with the potential off balance between the numbering of a pdb file and the amino acid sequence for the protein.
            #the try option with the continue option at the end save potential errors in case this trick is not necessary. 
            try:
                pos = int(float(raw_pos)) + offset
            except ValueError:
                continue

            #Translation between the amino acid code from different sources to align. 
            mut_3L = amino_acid_map.get(mut_1L, mut_1L)
            targets[pos]['muts'].add(mut_3L)
            targets[pos]['fallback_wt'] = amino_acid_map.get(orig_1L, "SER")

    print(f"Generating YAML files 1SNP con flexible structure(< {cutoff} Å) for{len(targets)} positions...")

    #Just setting up a count to numbering of files created. 
    count = 0
    
    #Thsi code selects to each sorted target position and amino acid the original wild type amino acid. 
    for pos, data in sorted(targets.items()):
        wt_aa = pdb_wt_map.get(pos, data['fallback_wt'])

        #This codes creates the protein configuration, adding the target residue to the list. 
        prot_configs = [format_residue_config(chain_id, pos, wt_aa, data['muts'])]

        #This code finds the neighbouring residues to be incorporated into the design yaml file as flexibles, according to the radious selected by user. 
        neighbors = find_3d_pocket_neighbors(pos, chain_id, pdb_data, cutoff_radius=cutoff)

        #A separate list is created to incorporate residues from the ligand or interacting protein. The mut_list=NONE allows for the setting of neighbours as flexible, without adding more information. 
        ligand_configs = []
        for n_chain, n_res, n_type in neighbors:
            if n_chain == chain_id:
                prot_configs.append(format_residue_config(n_chain, n_res, n_type, mut_list=None))
            else:
                ligand_configs.append(format_residue_config(n_chain, n_res, n_type, mut_list=None))

        #This code inserts the protein configuration for the mutations into the residue_configurations block of the yaml file. The join(prot_configs) puts together all the  blocks created corresponding to that file.  
        prot_block = "  residue_configurations:\n" + "".join(prot_configs)
        if "  residue_configurations: []" in clean_yaml:
            mutant_yaml = clean_yaml.replace("  residue_configurations: []", prot_block, 1)
        else:
            mutant_yaml = re.sub(r'protein:\s*\n', f'protein:\n{prot_block}', clean_yaml, count=1)

        #This code inserts the residues details for the ligand into the residue_configurations block of the ligand. 
        if ligand_configs:
            lig_block = "  residue_configurations:\n" + "".join(ligand_configs)
            mutant_yaml = re.sub(r'ligand:\s*\n\s*residue_configurations:\s*\[\]', f'ligand:\n{lig_block}', mutant_yaml)

        #Coding for naming of files created. 
        out_name = f"{yaml_prefix}.{chain_id}{pos}.yaml"
        with open(out_name, "w") as out_f:
            out_f.write(mutant_yaml)

        p_neigh = [f"{c}{r}" for c, r, _ in neighbors if c == chain_id]
        l_neigh = [f"{c}{r}" for c, r, _ in neighbors if c != chain_id]
        print(f"Generado {out_name} -> Diana: {chain_id}{pos} ({wt_aa}) | Shell: {len(p_neigh)} Protein, {len(l_neigh)} Ligand")
        count += 1

    print(f"\n{count} YAMLS files generated for 1SNP")


if __name__ == "__main__":
    print("--- 3D Flexible Mutation YAML Generator ---")
    
    # 1. Ask for the YAML file name
    yaml_input = input("Enter the name of your base YAML file (e.g., LptA-thanatin.yaml): ").strip()
    
    # 2. Ask for the radius
    radius_input = input("Enter the radius in Angstroms for flexible neighbors [Press Enter for 4.0]: ").strip()
    
    # 3. Safely handle the radius input
    if radius_input == "":
        cutoff_radius = 4.0
    else:
        try:
            cutoff_radius = float(radius_input)
        except ValueError:
            print("Invalid number entered. Defaulting to a 4.0 Å radius.")
            cutoff_radius = 4.0
            
    # 4. Run the generation process
    build_3d_flexible_1snp_yaml_files(
        base_yaml_path=yaml_input,
        cutoff=cutoff_radius
    )