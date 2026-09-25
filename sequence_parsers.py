"""

Sequence Parsers

================================================================================
THESIS TITLE : Analisis in silico del efecto de mutaciones en la interaccion de Thanatin con LptA 
AUTHOR       : Alexandro Mauricio Rizzo Gomez
MASTER       : Master in Bioinformatics
INSTITUTION  : Universidad Internacional de Valencia
SUPERVISORS  : Miquel Sendra
CREATED      : 2026-07-08
LAST UPDATE  : 2026-09-24
================================================================================
FILE PURPOSE:
    This code is part of a Master Thesis for a course in Bioinformatics at the International University of Valencia. 
    This code calculates the physicochemical differences between normal and mutated protein sequences by applying Grantham´s distance metrics.

================================================================================

"""

import pandas as pd
from Bio.Seq import Seq
#Imports the pandas and Seq modules for further processing. 

def extract_mutation_type(name):

#Parses the PoolParty row 'name' to identify the mutation class.

    if "1SNP" in name: 
        return "1SNP"
    if "2SNP" in name: 
        return "2SNP"
    if "del1bp" in name: 
        return "del1bp"
    if "del2bp" in name: 
        return "del2bp"
    if "del3bp" in name: 
        return "del3bp"
    if "ins_N_" in name or "ins1bp" in name: 
        return "ins1bp"
    if "ins_NN_" in name or "ins2bp" in name: 
        return "ins2bp"
    if "ins_NNN" in name or "ins3bp" in name: 
        return "ins3bp"
    return "Unknown"


def analyze_translation(dna_seq, wt_protein):
#This function translates DNA sequence and cleanly evaluates synonymous status and premature stop codons.
#accounts for length changes (deletions/insertions) so the natural stop codon is not misclassified as an early termination.

    clean_dna = dna_seq.replace("-", "")
    triplet_length = (len(clean_dna) // 3) * 3
    trimmed_dna = clean_dna[:triplet_length]
#First sequences are cleared of dash or empty strings, then the sequences is treated to leave only multuples of 3, to form sequences of only tripelts length for further analysis, the trimmed_dna variable.  
 
    mut_protein = str(Seq(trimmed_dna).translate(to_stop=False))
#Full translation including any internal stops deliveringa full amino acid sequence as a string.   
 
    wt_clean = wt_protein.rstrip("*")
#Removes the stop codon from teh very enf of the wild-type reference sequence for further analisis

    parts = mut_protein.split("*")
    protein_cleaned = parts[0]
# Split at the first stop codon encountered if there is one in the mutated sequence and keeps the amino acid sequence until that point for further analisis. 
  
    expected_len = len(wt_clean)
#This establsihes the baseline lenght of the normal protein, for further analisis.    
    
    if "*" in mut_protein:
#Check for stop codon in translation    
    
        first_stop_idx = mut_protein.index("*")
#Identifies the position of the stop codon.         
#If the stop happens substantially before the expected protein length, it's premature
#For del3bp, natural stop is at expected_len - 1
        has_new_stop = first_stop_idx < (expected_len - 2)
    else:
        has_new_stop = False
        
    is_synonymous = (protein_cleaned == wt_clean)
#Check if sequence is identical amino acid sequence to wild-type    
    
    return pd.Series([mut_protein, protein_cleaned, is_synonymous, has_new_stop])
#Stores the information into pandas series for integration into the dataframe. 


def parse_amino_acid_change(row, wt_protein):
  
#Parses exact 1-based absolute positions, wild-type amino acids, and mutant amino acids.
#Handles Point Mutations (1SNP/2SNP) and Indels (deletions, insertions, frameshifts).
    
    wt_prot = wt_protein
    mut_type = row["mutation_type"]
#Extracts reference protein and specific mutation classificacion from dataset
    
 #1. Handle Synonymous / Wild Type. If it is synonumours, stores none value for position, original and mutant columns. 
    if row["is_synonymous"]:
        return pd.Series(["None", "None", "None"])
    
#2. Handle Point Mutations (1SNP & 2SNP): compare full-length translation
    if mut_type in ["1SNP", "2SNP"]:
        clean_dna = row["seq"].replace("-", "")
        triplet_length = (len(clean_dna) // 3) * 3
        full_mut_prot = str(Seq(clean_dna[:triplet_length]).translate(to_stop=False))
        
        pos_list, wt_list, mut_list = [], [], []
#Creates lists for temporary storage for the mutation details(position, original amino acid and mutant amino acid).
        
        for i, (orig, mut) in enumerate(zip(wt_prot, full_mut_prot), start=1):
            if orig != mut:
                pos_list.append(str(i))
                wt_list.append(orig)
                mut_list.append(mut if mut != "*" else "* (Stop)")
# Index i starts at 1, corresponding to absolute protein residue number. Through the enumerate and zip functions aligns 
#the wild type and mutant sequence and extracts information of mutations.                 
                
        return pd.Series([",".join(pos_list), ",".join(wt_list), ",".join(mut_list)])
#The temporary lists are included as contests into the dataset.
    
# 3. Handle Indels (del1bp, del2bp, del3bp, ins1bp, ins2bp, ins3bp, etc.)
#Pulls truncated mutant protein sequences. 
    mut_prot = row["mut_protein_clean"]
    first_diff = 1
#Tracks position at which sequences begin to differ. 
    min_len = min(len(wt_prot), len(mut_prot))    
    for i in range(min_len):
#Ensuring iteration occurs through the shorters sequence lenght. 
        
        if wt_prot[i] != mut_prot[i]:
            first_diff = i + 1
            break
#If it finds mitchmatches in amino acid sequence between wildtype and mutant, tracks the position.

    else:
        first_diff = min_len + 1

    len_diff = len(mut_prot) - len(wt_prot)
#Calculates lenght of mutated sequence to determine change in amino acids count.     
    
    if len_diff < 0:  
#Selection for frameshift Deletion

        del_len = abs(len_diff)
#Calculation of missing amino acids. 
        wt_sub = wt_prot[first_diff - 1 : first_diff - 1 + del_len]
        mut_sub = f"Δ{del_len}AA"
        pos_str = f"{first_diff}-{first_diff - 1 + del_len}" if del_len > 1 else f"{first_diff}"
        return pd.Series([pos_str, wt_sub, mut_sub])
#Code for determine details of the mutation   
        
    elif len_diff > 0:  
#Selection for frameshift Insertion

        ins_len = len_diff
#Calculation for amino acid gain
        wt_sub = "-"
        mut_sub = mut_prot[first_diff - 1 : first_diff - 1 + ins_len]
        pos_str = f"ins_{first_diff}"
        return pd.Series([pos_str, wt_sub, mut_sub])
#Code for determine details of the mutation 
 
    else:  
# Selection for frameshift without overall length change (or leading to truncation)

        wt_sub = wt_prot[first_diff - 1 : first_diff] if first_diff <= len(wt_prot) else "-"
        mut_sub = "FS_Truncation"
        return pd.Series([f"{first_diff}_frameshift", wt_sub, mut_sub])
#Code for determine details of the mutation

def parse_double_snp_split(row, wt_protein):
   
#Parses a 2SNP variant into a dictionary with keys matching 
#pos1, orig1, mut1, pos2, orig2, mut2.
    
    wt_prot = wt_protein
    clean_dna = row["seq"].replace("-", "")
    triplet_length = (len(clean_dna) // 3) * 3
    full_mut_prot = str(Seq(clean_dna[:triplet_length]).translate(to_stop=False))
#Preprocessing of the sequences and translation into amino acid sequence. 
 
    positions, origs, muts = [], [], []
#Creates lists for temporary storage for the mutation details(position, original amino acid and mutant amino acid).
    
    min_len = min(len(wt_prot), len(full_mut_prot))
#Ensuring iteration occurs through the shorters sequence lenght.     
    
    for i in range(min_len):
        orig = wt_prot[i]
        mut = full_mut_prot[i]
        if orig != mut:
            positions.append(i + 1)  
            origs.append(orig)
            muts.append(mut if mut != "*" else "*")
#Code for determine details of the mutation          
         
#Always return a dictionary with the 6 expected keys, equal to two mutations. 
    if len(positions) >= 2:
        return {
            "pos1": positions[0], "orig1": origs[0], "mut1": muts[0],
            "pos2": positions[1], "orig2": origs[1], "mut2": muts[1]
        }
    else:
        return {
            "pos1": None, "orig1": None, "mut1": None,
            "pos2": None, "orig2": None, "mut2": None
        }
#Code to treat potential problems with the sequence and avoid further problems down the line byt returning a None value. 