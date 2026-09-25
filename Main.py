"""

Main.py

================================================================================
THESIS TITLE : Analisis in silico del efecto de mutaciones en la interaccion de Thanatin con LptA 
AUTHOR       : Alexandro Mauricio Rizzo Gomez
MASTER       : Master in Bioinformatics
INSTITUTION  : Universidad Internacional de Valencia
SUPERVISORS  : Miquel Sendra
CREATED      : 2026-06-28
LAST UPDATE  : 2026-09-24
================================================================================
FILE PURPOSE:
    This code is part of a Master Thesis for a course in Bioinformatics at the International University of Valencia. 
    This code is the main body of the Mutantgenerator python program of the thesis. It handles the user interaction, 
    coordinating the generation of mutant DNA sequences with the support of the PoolParty library, processes the 
    translation, score the mutations and exports the final datasets and reports. 

================================================================================

"""
#Imports the required libraries for the correct functioning of the program
import itertools
import math
import os
import pandas as pd
import poolparty as pp
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Seq import Seq

#These functions are imported from other pieces of code generated in this workflow. 
from utils.grantham import calculate_grantham_score
from utils.probability import classify_nuc_mutation, calculate_ts_tv_probability
from utils.sequence_parsers import (
    analyze_translation,
    extract_mutation_type,
    parse_amino_acid_change,
    parse_double_snp_split,
)


def get_user_inputs():
    print("=" * 60)
    print("   MUTANTGENERATOR, IN SILICO MUTAGENESIS")
    print("=" * 60)

    #1. Project Identifier. The user is asked to name the mutagenesis run for identification. 
    project_name = input("Enter project/run name [e.g., BamA_Darobactin]: ").strip()
    if not project_name:
        project_name = "Mutagenesis_Project"

    #2. Complete Open Reading Frame (ORF) with multiline paste support, The code asks user for the original genomic sequence to mutate. 
    print("\nPaste your nucleotide sequence below (press Enter twice on an empty line when finished):")
    lines = []
    #Starts the empy list that stores the individual lines of input. There is a safety mechanism, to avoid crashes, where the system detects errors 
    #due to large blocks of text or wrong keyboards pressed. The code requires double pressing Enter to confirm end of pasting content, only if lines
    #have been filled. The .strip function deletes white spaces. 
    while True:
        try:
            line = input()
            if not line.strip():
                if lines:
                    break
                continue
            lines.append(line.strip())
        except EOFError:
            break

    raw_seq = "".join(lines).upper()
    clean_seq = "".join([c for c in raw_seq if c in "ATGCN"])
    print(f"\nLoaded sequence of length: {len(clean_seq)} bp")

    #3. Targeted Regions. Prompts the user to determine how many regions are to be mutagenized and the start and end positions of those regions. 
    while True:
        try:
            num_regions = int(input("\nHow many distinct regions do you want to target? "))
            if num_regions >= 1:
                break
            print("Enter an integer >= 1.")
        except ValueError:
            print("Invalid input. Please enter a valid number.")

    regions = []
    print(f"\nDefine boundaries for {num_regions} target region(s).")
    print("Note: Input 1-based nucleotide positions (e.g., 88 to 135).")
    for r in range(1, num_regions + 1):
        while True:
            try:
                start_pos = int(input(f" - Region {r} START nucleotide (1-based): "))
                end_pos = int(input(f" - Region {r} END nucleotide (1-based): "))
                if 1 <= start_pos < end_pos <= len(clean_seq):
                    # Convert to 0-indexed slicing coordinates
                    regions.append((start_pos - 1, end_pos, f"R{r}"))
                    break
                print(f"Boundaries must satisfy 1 <= START < END <= {len(clean_seq)}.")
            except ValueError:
                print("Invalid integer input. Please try again.")

    # 4. Mutation Selection of Top-N Export. Use is asked what kind of mutations does he want to perform and how many sequences does he want the program to extract for further analysis.  
    print("\nSelect mutation types and specify candidate export:")
    
    #1SNP selection & quantity selection
    gen_1snp = input(" - Include single point mutations (1SNP)? [Y/n]: ").strip().lower() != "n"
    top_1snp = 50
    if gen_1snp:
        val = input("   How many top 1SNP candidates to export? [default: 50]: ").strip()
        top_1snp = int(val) if val.isdigit() and int(val) > 0 else 50

    #2SNP selection & quantity selection
    gen_2snp = input(" - Include double point mutations (2SNP)? [Y/n]: ").strip().lower() != "n"
    top_2snp = 50
    if gen_2snp:
        val = input("   How many top 2SNP candidates to export? [default: 50]: ").strip()
        top_2snp = int(val) if val.isdigit() and int(val) > 0 else 50
    
    #InDel selection & quantity selection
    gen_indels = input(" - Include InDels (insertions/deletions)? [y/N]: ").strip().lower() == "y"
    indel_sizes = []
    top_indels = 50
    if gen_indels:
        sizes_in = input("   Enter InDel sizes to scan separated by commas [e.g., 1,2,3]: ").strip()
        for size in sizes_in.split(","):
            size = size.strip()
            if size in ["1", "2", "3"]:
                indel_sizes.append(int(size))
        if not indel_sizes:
            indel_sizes = [1, 2, 3]
        val = input("   How many top candidates to export per InDel sub-type? [default: 50]: ").strip()
        top_indels = int(val) if val.isdigit() and int(val) > 0 else 50

    #5. Experimental Modeling Mode. User can determine the type of mutation context to simulate. 
    print("\nSelect Ts/Tv Mutation Model background:")
    print(" [1] genomic_wildtype (basal rate 1e-9)")
    print(" [2] plasmid_epPCR (error-prone PCR rate 1e-3)")
    mode_choice = input("Enter choice [1/2, default 1]: ").strip()
    exp_mode = "plasmid_epPCR" if mode_choice == "2" else "genomic_wildtype"

    #Program returns input as a dictionary. 
    return {
        "project_name": project_name,
        "orf": clean_seq,
        "regions": regions,
        "gen_1snp": gen_1snp,
        "top_1snp": top_1snp,
        "gen_2snp": gen_2snp,
        "top_2snp": top_2snp,
        "gen_indels": gen_indels,
        "indel_sizes": indel_sizes,
        "top_indels": top_indels,
        "exp_mode": exp_mode,
    }


def main():
    #The main function of the program runs the get_user_inputs function. 
    cfg = get_user_inputs()
    orf = cfg["orf"]
    project_name = cfg["project_name"]
    regions = cfg["regions"]

    # Analysis of the wildtype sequence, ensuring the length is ok for further analysis as well as determining variables for later use. 
    triplet_len = (len(orf) // 3) * 3
    wt_protein = str(Seq(orf[:triplet_len]).translate(to_stop=False)).rstrip("*")
    total_target_len = sum(end - start for start, end, _ in regions)
    total_context_len = len(orf)

    print(f"\n Initializing library generation for {project_name}...")
    print(f"Gene length: {len(orf)} bp | Protein baseline: {len(wt_protein)} aa")

    # -------------------------------------------------------------------
    # 1. GENERATE POOLPARTY LIBRARIES PER REGION
    # -------------------------------------------------------------------
    region_dfs = []

    #For each region selected by user, the code cleave steh target sequence for further analysis with Poolparty
    for start_idx, end_idx, reg_label in regions:
        prefix_seq = orf[:start_idx]
        target_seq = orf[start_idx:end_idx]
        suffix_seq = orf[end_idx:]

        # Creates an empty list to act as a temporary storage.  
        snp_dfs = []

        #If use has selected mutations of 1SNP, a new fresh PoolParty environment is initiated, target_seq is feed into the 
        #mutagenize method and with the hard coded num_mutationsm, sequential mode, it generates a library with nym_cycles= 1 and seed=42, 
        #storing the DNA library as panda data frames df_1. It does a similar job with the 2NSP mutagenesis if the user has request it. 
        if cfg["gen_1snp"]:
            pp.init()
            df_1 = pp.from_seq(target_seq).mutagenize(
                num_mutations=1, prefix=f"1SNP_{reg_label}", mode="sequential"
            ).generate_library(num_cycles=1, seed=42)
            snp_dfs.append(df_1)

        if cfg["gen_2snp"]:
            pp.init()
            df_2 = pp.from_seq(target_seq).mutagenize(
                num_mutations=2, prefix=f"2SNP_{reg_label}", mode="sequential"
            ).generate_library(num_cycles=1, seed=42)
            snp_dfs.append(df_2)

        if snp_dfs:
            df_snps = pd.concat(snp_dfs, ignore_index=True)
        else:
            df_snps = pd.DataFrame(columns=["name", "seq"])

        # InDels generation. Creating the indel_dfs for temporary storage, and with the variable size according to the user input,the code first 
        #initiates a new poolparty section, and with the deletion_scan function from pool party, generates the library/ies according to user input,
        #going through the size of deleletion or insertions, providing each squence the descriptive  prefix. Insertions are carried out in a similar fashion.
        #All sequences are added to the temporary list. 
        indel_dfs = []
        if cfg["gen_indels"]:
            for size in cfg["indel_sizes"]:
                pp.init()
                d_scan = pp.from_seq(target_seq).deletion_scan(
                    deletion_length=size, prefix=f"del{size}bp_{reg_label}", mode="sequential"
                ).generate_library(num_cycles=1, seed=42)
                indel_dfs.append(d_scan)

                pp.init()
                ins_motif = "N" * size
                i_scan = pp.from_seq(target_seq).insertion_scan(
                    insertion_pool=pp.from_seq(ins_motif), prefix=f"ins{size}bp_{reg_label}", mode="sequential"
                ).generate_library(num_cycles=1, seed=42)
                indel_dfs.append(i_scan)

        #All of the indel mutations are put together with the dataframe df_snps, that has the data of all snps, into a unified data table, reg_combined, 
        #with the concat function from pandas.
        reg_combined = pd.concat([df_snps] + indel_dfs, ignore_index=True)
        
        #The next piece of code re-stitch mutated fragment back into full ORF context. With the apply function from pandas, lambda assigns to every row of seq the variable frag, 
        #cleans it of any dash characters, and inserts the previos prefix and suffix to the fragments to create the full mutated sequence that is incorporated in the regions_dfs 
        #dataframe with the append function 
        reg_combined["seq"] = reg_combined["seq"].apply(
            lambda frag: prefix_seq + frag.replace("-", "") + suffix_seq
        )
        region_dfs.append(reg_combined)

    #Finally, all the different lists are put together into a dataframe that incorporates all the information. It also has a quick check in case mutations are not being generated 
    #through a mistake in the system, sending a error if the dataframe is empty.     
    df = pd.concat(region_dfs, ignore_index=True)
    if df.empty:
        print("No mutations generated based on current configuration. Exiting.")
        return

    # -------------------------------------------------------------------
    # 2. TRANSLATION & REDUNDANCY ANALYSIS (FULL DNA LIBRARY LEVEL)
    # -------------------------------------------------------------------
    #This section deals with the translation of the mutated sequences into amino acid sequences, grouping them according to similitude. 
    print("Translating sequences and running biological filters...")
    
    #With the function apply every styring of seq is parsed with the created analyze_translation creating the information to fill the dataframe, with the 
    #mutated protein sequence and its characteristics regarding stop codon or synonymous status.
    
    df[["mut_protein", "mut_protein_clean", "is_synonymous", "has_stop_codon"]] = df["seq"].apply(
        lambda dna_string: analyze_translation(dna_string, wt_protein)
    )

    #This piece of coding extraxcts the type of mutations from the name on the tags with the extract_mutation_type function. 
    df["mutation_type"] = df["name"].apply(extract_mutation_type)
    #This sections compares the difference between the sequence to clasify the type of nucleotide mutation for further analysis for translations transversions. 
    df["nuc_mutation_detail"] = df.apply(lambda r: classify_nuc_mutation(r, orf), axis=1)
    #The next lines identify the type of mutation, extracting the information of the original amino acid, the mutant one and the position. 
    df[["aa_position", "wt_amino_acid", "mut_amino_acid"]] = df.apply(
        lambda row: parse_amino_acid_change(row, wt_protein), axis=1
    )

    #This piece of code tracks the redundancy across all generated DNA secuence, counting the unique secuences and storign the values with the to_dict function. With the map funtion 
    #appied to protein_counts to assign the equivalent_sequences_count. 
    protein_counts = df["mut_protein_clean"].value_counts().to_dict()
    df["equivalent_sequences_count"] = df["mut_protein_clean"].map(protein_counts)
    #This code assigns ID to each unique protein, through the enumarate and unique functions, with a formating for the index  and with the map function assiging to each identical 
    #protein sequence the same generated protein group ID GRP_...
    unique_proteins = {p: f"GRP_{i+1:05d}" for i, p in enumerate(df["mut_protein_clean"].unique())}
    df["protein_group_id"] = df["mut_protein_clean"].map(unique_proteins)

    # -------------------------------------------------------------------
    # 3. MULTI-PARAMETRIC SCORING
    # -------------------------------------------------------------------
    #In this section, the scoring of the amino acid sequence according to the BLOSUM62, Grantham distance and probability is generated.
    print("Calculating BLOSUM62, Grantham, and Ts/Tv probabilities...")
    
    #Initialization of the PairwiseAligner from Biopython library
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    #Code for the calculation of the BLOSUM62 score, with the apply method between the wildtype protein and the mutant one (m). The if m else 0.0 is to prevent errors due to empy 
    #strings from early stop codons by avoiding the parsing of empty sequences.  
    df["blosum62_scores"] = df["mut_protein_clean"].apply(
        lambda m: aligner.score(wt_protein, m) if m else 0.0
    )
    
    #Calculation of Grantham distance with the function coded in the separate file. 
    df["grantham_scores"] = df["mut_protein_clean"].apply(
        lambda m: calculate_grantham_score(wt_protein, m)
    )
    
    #Runs the mutation probability calculation. Detail corresponds to the detailed nucleotide mutation type,i.e. A to C or A to T...  
    df["mutation_probability"] = df["nuc_mutation_detail"].apply(
        lambda m: calculate_ts_tv_probability(m, total_target_len, total_context_len, cfg["exp_mode"])
    )
    #Generation of the amino acid sequence column, prior to future exports coding.  
    df["amino_acid_sequence"] = df["mut_protein_clean"]

    # -------------------------------------------------------------------
    # 4. EXPORT MASTER DATASET (PRESERVING COMPLETE LIBRARY REDUNDANCY)
    # -------------------------------------------------------------------
    #This section produces the production of the csv master file generated. 
    
    #Code to the generation of the file name, base on project name by user, and exporting the correct columns into it. 
    master_csv = f"{project_name}_mutatgenesis.csv"
    export_columns = [
        "protein_group_id", "equivalent_sequences_count",
        "aa_position", "wt_amino_acid", "mut_amino_acid", "mutation_type",
        "blosum62_scores", "grantham_scores", "mutation_probability", "amino_acid_sequence"
    ]
    df[export_columns].to_csv(master_csv, index=False)
    print(f"Exported complete master dataset ({len(df):,} rows) to '{master_csv}'")

    # -------------------------------------------------------------------
    # 5. DEDUPLICATION & DYNAMIC TOP-N EXPORTS PER SUB-TYPE
    # -------------------------------------------------------------------
    #This section organized the creation and export of sequences selected by the user.
    
    #Filtration of synonymous and broken sequences of amino acids and copying them into a new dataframe
    clean_variants = df[(~df["is_synonymous"]) & (~df["has_stop_codon"])].copy()

    #Sorting of the sequences according to mutation probability and elimination of duplications for further analysis. 
    clean_dedup = clean_variants.sort_values(by="mutation_probability", ascending=False)
    clean_dedup = clean_dedup.drop_duplicates(subset=["mut_protein_clean"], keep="first").copy()
    print(f"Deduplicated viable variants: {len(clean_variants):,} DNA variants -> {len(clean_dedup):,} unique protein outcomes.")

    #Handling of top 1SNP mutants. If use has selected a specific number, or otherwise with the standard 50 option, the code
    #sorts and selects the top disruptive mutations according to Grantham distances and BLOSUM, and exporting key parameters for further analisis.
    if cfg["gen_1snp"]:
        n_1snp = cfg["top_1snp"]
        df_1snp = clean_dedup[clean_dedup["mutation_type"] == "1SNP"].sort_values(
            by=["grantham_scores", "blosum62_scores"], ascending=[False, True]
        ).head(n_1snp).copy()
        
        df_1snp_export = df_1snp[["aa_position", "wt_amino_acid", "mut_amino_acid"]].rename(
            columns={"aa_position": "position", "wt_amino_acid": "original", "mut_amino_acid": "mutant"}
        )
        df_1snp_export.to_csv("1snp_mutants.csv", index=False)
        print(f"Exported top {len(df_1snp_export)} unique 1SNP candidates to '1snp_mutants.csv'")

    #Handling of top 2SNP mutants. 
    if cfg["gen_2snp"]:
        n_2snp = cfg["top_2snp"]
        df_2snp = clean_dedup[clean_dedup["mutation_type"] == "2SNP"].sort_values(
            by=["grantham_scores", "blosum62_scores"], ascending=[False, True]
        ).head(n_2snp).copy()
        
        
        #After a basic check to prevent errors due to empty data, the code produces with the parse_double_snp_split a list of double changes in the 
        #amino acid sequences storing the changes in two different occassions for further processing, droping intracodon mutations from the final result.
        #.astype(int) is a mechanisms to revert the format to the position columns when dropping lines with empty columns.          
        if not df_2snp.empty:
            parsed_records = [parse_double_snp_split(row, wt_protein) for _, row in df_2snp.iterrows()]
            df_2snp_export = pd.DataFrame(parsed_records).dropna().copy()
            if not df_2snp_export.empty:
                df_2snp_export["pos1"] = df_2snp_export["pos1"].astype(int)
                df_2snp_export["pos2"] = df_2snp_export["pos2"].astype(int)
                df_2snp_export.to_csv("2snp_mutants.csv", index=False)
                print(f"Exported top {len(df_2snp_export)} unique 2SNP candidates to '2snp_mutants.csv'")

    #Top InDel per sub-type, if user has determined this options. By filtering by indel type and size, the code sorts them and creates the file to export. 
    #It also manages the potential that the mutations does not generate valid sequences. 
    if cfg["gen_indels"]:
        n_indels = cfg["top_indels"]
        indel_export_cols = [
            "aa_position", "wt_amino_acid", "mut_amino_acid", 
            "mutation_type", "grantham_scores", "blosum62_scores", 
            "mutation_probability", "amino_acid_sequence"
        ]

        for size in cfg["indel_sizes"]:
            # Deletions
            del_type = f"del{size}bp"
            df_del_sub = clean_dedup[clean_dedup["mutation_type"] == del_type].sort_values(
                by=["grantham_scores", "blosum62_scores"], ascending=[False, True]
            ).head(n_indels).copy()

            del_file = f"{del_type}_mutants.csv"
            if not df_del_sub.empty:
                df_del_sub[indel_export_cols].to_csv(del_file, index=False)
                print(f"Exported top {len(df_del_sub)} unique {del_type} candidates to '{del_file}'")
            else:
                print(f"No viable/non-stop candidates found for {del_type}.")

            # Insertions
            ins_type = f"ins{size}bp"
            df_ins_sub = clean_dedup[clean_dedup["mutation_type"] == ins_type].sort_values(
                by=["grantham_scores", "blosum62_scores"], ascending=[False, True]
            ).head(n_indels).copy()

            ins_file = f"{ins_type}_mutants.csv"
            if not df_ins_sub.empty:
                df_ins_sub[indel_export_cols].to_csv(ins_file, index=False)
                print(f"Exported top {len(df_ins_sub)} unique {ins_type} candidates to '{ins_file}'")
            else:
                print(f"No viable/non-stop candidates found for {ins_type}.")

    # -------------------------------------------------------------------
    # 6. NUMERICAL EXPORTS FOR GRAPHING & DOWNSTREAM STATISTICS
    # -------------------------------------------------------------------
    print(f"Generating numerical dataframes for downstream plotting...")

    # A. Descriptive Statistics by Mutation Type (Grantham, BLOSUM62, Ts/Tv)
    metrics = ["grantham_scores", "blosum62_scores", "mutation_probability"]
    stat_records = []

    #By selecting type of mutation, and for each variable selected in the metrics, the describe function from pandas calculates the values, 
    #which are then saved into a dictionary through the append function and eventually added to the statistical summary file.  
    for m_type, group in clean_dedup.groupby("mutation_type"):
        for metric in metrics:
            desc = group[metric].describe()
            stat_records.append({
                "mutation_type": m_type,
                "metric": metric,
                "count": desc["count"],
                "mean": desc["mean"],
                "std": desc["std"],
                "min": desc["min"],
                "25%": desc["25%"],
                "median": desc["50%"],
                "75%": desc["75%"],
                "max": desc["max"]
            })

    stats_df = pd.DataFrame(stat_records)
    stats_csv = f"{project_name}_statistical_summary.csv"
    stats_df.to_csv(stats_csv, index=False)
    print(f"Exported numerical distribution stats to '{stats_csv}'")

    # B. Amino Acid Mutation Frequencies (for substitution matrices & heatmaps)
    #This piece of code groups mutations by the specific type of mutation of amino acid, counting (size option) and convert it that count into a columnd through the reset_index
    aa_freq_df = (
        clean_dedup.groupby(["mutation_type", "wt_amino_acid", "mut_amino_acid"])
        .size()
        .reset_index(name="count")
    )
    
    #The values are then transformed into frequencies for further analyis and exproted to a new file. 
    type_totals = aa_freq_df.groupby("mutation_type")["count"].transform("sum")
    aa_freq_df["relative_frequency"] = aa_freq_df["count"] / type_totals
    aa_freq_df = aa_freq_df.sort_values(by=["mutation_type", "count"], ascending=[True, False])

    freq_csv = f"{project_name}_aa_mutation_frequencies.csv"
    aa_freq_df.to_csv(freq_csv, index=False)
    print(f"Exported amino acid substitution frequencies to '{freq_csv}'")

    # -------------------------------------------------------------------
    # 7. GENERATE COMPREHENSIVE TEXT SUMMARY REPORT
    # -------------------------------------------------------------------
    # Final section that generates summary report. 
    summary_file = f"{project_name}_generation_summary.txt"
    print(f"Writing comprehensive summary to '{summary_file}'...")

    df_valid_1snp = clean_dedup[clean_dedup["mutation_type"] == "1SNP"]
    df_valid_2snp = clean_dedup[clean_dedup["mutation_type"] == "2SNP"]

    with open(summary_file, "w") as f:
        f.write("=================================================\n")
        f.write(f"      {project_name} Mutational Library Report    \n")
        f.write("=================================================\n\n")
        
        # Section 1: Overall Generation Metrics
        f.write("--- 1. OVERALL GENERATION METRICS ---\n")
        f.write(f"Total DNA sequences generated: {len(df):,}\n")
        f.write(f"Sequences with early stop codons (Nonsense): {df['has_stop_codon'].sum():,}\n")
        f.write(f"Synonymous sequences (Silent / Equivalent to WT): {df['is_synonymous'].sum():,}\n")
        f.write(f"Viable Missense sequences (Total DNA variants): {len(clean_variants):,}\n")
        f.write(f"Viable Missense sequences (Unique Protein variants): {len(clean_dedup):,}\n")
        f.write(f"Total unique translated protein variants: {df['mut_protein_clean'].nunique():,}\n\n")

        # Section 2: Breakdown by Mutation Type (Total)
        f.write("--- 2. BREAKDOWN BY MUTATION TYPE (Total) ---\n")
        type_counts = df['mutation_type'].value_counts()
        for m_type, count in type_counts.items():
            f.write(f" - {m_type:8s}: {count:,}\n")
        f.write("\n")

        # Section 3: Sequence Equivalence & Redundancy Analysis
        f.write("--- 3. SEQUENCE EQUIVALENCE & REDUNDANCY ANALYSIS ---\n")
        redundant_groups = df[df["equivalent_sequences_count"] > 1]
        f.write(f"Number of DNA sequences that share an identical protein outcome: {len(redundant_groups):,}\n")
        f.write(f"Number of distinct redundant protein groups: {redundant_groups['mut_protein_clean'].nunique():,}\n\n")
        
        f.write("Top 10 Most Redundant Protein Outcomes:\n")
        top_redundant = (
            df.groupby(["protein_group_id", "aa_position", "wt_amino_acid", "mut_amino_acid"])
            .size()
            .reset_index(name="count")
            .sort_values(by="count", ascending=False)
            .head(10)
        )
        for _, r in top_redundant.iterrows():
            f.write(f" - Group {r['protein_group_id']} ({r['wt_amino_acid']}{r['aa_position']}{r['mut_amino_acid']}): generated by {r['count']} distinct DNA sequences\n")
        f.write("\n")

        # Section 4: 1SNP Specific Substitution Tracking
        f.write("--- 4. 1SNP AMINO ACID SUBSTITUTION TRACKING ---\n")
        f.write("(Showing valid unique missense 1SNPs across all targeted regions)\n\n")
        
        aa_changes_1snp = df_valid_1snp.groupby(['wt_amino_acid', 'mut_amino_acid']).size().reset_index(name='count')
        aa_changes_1snp = aa_changes_1snp.sort_values(by='count', ascending=False)
        
        if not aa_changes_1snp.empty:
            for _, row in aa_changes_1snp.iterrows():
                f.write(f"{row['wt_amino_acid']} --> {row['mut_amino_acid']} : {row['count']} generated\n")
        else:
            f.write("No valid 1SNP substitutions found.\n")
        f.write("\n")

        # Section 5: 2SNP Specific Substitution Tracking
        f.write("--- 5. 2SNP AMINO ACID SUBSTITUTION TRACKING ---\n")
        f.write("(Showing valid unique missense 2SNPs across all targeted regions)\n\n")
        
        aa_changes_2snp = df_valid_2snp.groupby(['wt_amino_acid', 'mut_amino_acid']).size().reset_index(name='count')
        aa_changes_2snp = aa_changes_2snp.sort_values(by='count', ascending=False)
        
        if not aa_changes_2snp.empty:
            for _, row in aa_changes_2snp.iterrows():
                f.write(f"[{row['wt_amino_acid']}] --> [{row['mut_amino_acid']}] : {row['count']} generated\n")
        else:
            f.write("No valid 2SNP substitutions found.\n")
        f.write("\n")

        # Section 6: Biochemical & Evolutionary Scoring (Descriptive Statistics)
        f.write("--- 6. BIOCHEMICAL & EVOLUTIONARY SCORING (DESCRIPTIVE STATISTICS) ---\n")
        
        if not df_valid_1snp.empty:
            f.write("\n>> 1SNP STATISTICS:\n")
            gran_desc = df_valid_1snp["grantham_scores"].describe()
            f.write("  Grantham Scores (Biochemical Distance):\n")
            f.write(f"    Mean  : {gran_desc['mean']:.2f}\n")
            f.write(f"    Median: {gran_desc['50%']:.2f}\n")
            f.write(f"    StdDev: {gran_desc['std']:.2f}\n")
            f.write(f"    Min   : {gran_desc['min']:.0f}\n")
            f.write(f"    Max   : {gran_desc['max']:.0f}\n")
            if gran_desc["mean"] > 100:
                f.write("    -> Interpretation: The 1SNP library leans towards RADICAL biochemical changes on average.\n")
            else:
                f.write("    -> Interpretation: The 1SNP library leans towards CONSERVATIVE biochemical changes on average.\n")
                
            blo_desc = df_valid_1snp["blosum62_scores"].describe()
            f.write("\n  BLOSUM62 Scores (Evolutionary Viability):\n")
            f.write(f"    Mean  : {blo_desc['mean']:.2f}\n")
            f.write(f"    Median: {blo_desc['50%']:.2f}\n")
            f.write(f"    StdDev: {blo_desc['std']:.2f}\n")
            f.write(f"    Min   : {blo_desc['min']:.0f}\n")
            f.write(f"    Max   : {blo_desc['max']:.0f}\n")
            
            prob_desc = df_valid_1snp["mutation_probability"].describe()
            f.write("\n  Ts/Tv Mutation Probability (Context Background):\n")
            f.write(f"    Mean  : {prob_desc['mean']:.4e}\n")
            f.write(f"    Median: {prob_desc['50%']:.4e}\n")
            f.write(f"    StdDev: {prob_desc['std']:.4e}\n")
            f.write(f"    Min   : {prob_desc['min']:.4e}\n")
            f.write(f"    Max   : {prob_desc['max']:.4e}\n")

        if not df_valid_2snp.empty:
            f.write("\n-------------------------------------------------\n")
            f.write(">> 2SNP STATISTICS:\n")
            gran_desc2 = df_valid_2snp["grantham_scores"].describe()
            f.write("  Grantham Scores (Biochemical Distance):\n")
            f.write(f"    Mean  : {gran_desc2['mean']:.2f}\n")
            f.write(f"    Median: {gran_desc2['50%']:.2f}\n")
            f.write(f"    StdDev: {gran_desc2['std']:.2f}\n")
            f.write(f"    Min   : {gran_desc2['min']:.0f}\n")
            f.write(f"    Max   : {gran_desc2['max']:.0f}\n")
            if gran_desc2["mean"] > 100:
                f.write("    -> Interpretation: The 2SNP library leans towards RADICAL biochemical changes on average.\n")
            else:
                f.write("    -> Interpretation: The 2SNP library leans towards CONSERVATIVE biochemical changes on average.\n")
                
            blo_desc2 = df_valid_2snp["blosum62_scores"].describe()
            f.write("\n  BLOSUM62 Scores (Evolutionary Viability):\n")
            f.write(f"    Mean  : {blo_desc2['mean']:.2f}\n")
            f.write(f"    Median: {blo_desc2['50%']:.2f}\n")
            f.write(f"    StdDev: {blo_desc2['std']:.2f}\n")
            f.write(f"    Min   : {blo_desc2['min']:.0f}\n")
            f.write(f"    Max   : {blo_desc2['max']:.0f}\n")
            
            prob_desc2 = df_valid_2snp["mutation_probability"].describe()
            f.write("\n  Ts/Tv Mutation Probability (Context Background):\n")
            f.write(f"    Mean  : {prob_desc2['mean']:.4e}\n")
            f.write(f"    Median: {prob_desc2['50%']:.4e}\n")
            f.write(f"    StdDev: {prob_desc2['std']:.4e}\n")
            f.write(f"    Min   : {prob_desc2['min']:.4e}\n")
            f.write(f"    Max   : {prob_desc2['max']:.4e}\n\n")

        # Section 7: Top Mutated Residues
        f.write("--- 7. TOP 5 MOST MUTATED POSITIONS (1SNP) ---\n")
        if not df_valid_1snp.empty:
            top_positions = df_valid_1snp["aa_position"].value_counts().head(5)
            for pos, count in top_positions.items():
                f.write(f" - Position {pos}: {count} unique valid mutations\n")
        f.write("\n")

        # Section 8: InDel Summary
        if cfg["gen_indels"]:
            f.write("--- 8. INDEL SUMMARY (1-3 bp Insertions and Deletions) ---\n")
            indel_df = df[df["mutation_type"].str.contains("del|ins", regex=True)]
            if not indel_df.empty:
                summary_table = indel_df.groupby(["mutation_type", "has_stop_codon"]).size().unstack(fill_value=0)
                summary_table.columns = ["Non-Stop", "Premature Stop"]
                for m_type, row_vals in summary_table.iterrows():
                    f.write(f" - {m_type:8s} -> In-Frame/Viable: {row_vals['Non-Stop']:,} | Truncated: {row_vals['Premature Stop']:,}\n")

    print(f"Full analysis and summary report exported to '{summary_file}'")


if __name__ == "__main__":
    main()