"""
Transition/Transversion (Ts/Tv) Classification & Binomial Mutation Models for E. coli.

================================================================================
THESIS TITLE : Analisis in silico del efecto de mutaciones en la interaccion de Thanatin con LptA 
AUTHOR       : Alexandro Mauricio Rizzo Gomez
MASTER       : Master in Bioinformatics
INSTITUTION  : Universidad Internacional de Valencia
SUPERVISORS  : Miquel Sendra
CREATED      : 2026-07-10
LAST UPDATE  : 2026-09-24
================================================================================
FILE PURPOSE:
    This code is part of a Master Thesis for a course in Bioinformatics at the International University of Valencia. 
    This code calculates the probability of a type of mutation occurring according to a binomial combination model.

================================================================================

"""
import math
# Import math module for calculation of binomial combinations

TRANSITIONS = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
#Definition of transitions (purine to purine, pyrimide to pyrimidine)for later analysis. If a mutation is not present, it will be considered a transversion(i.e. "A" "T").

def classify_nuc_mutation(row, wt_dna_seq):
# Funcation that analyses type of mutation that has occurred.    
    wt_seq = wt_dna_seq
    mut_seq = row["seq"].replace("-", "")#elimination of any dash characters present in secuences as result of operations performed on it
    
    if len(wt_seq) != len(mut_seq):
        return "Indel"
#Check for indels according to sequence lenght, if it differs, it is a result of a indel   
 
    changes = []
    for orig, mut in zip(wt_seq, mut_seq):
        if orig != mut:
            sub_type = "Ts" if (orig, mut) in TRANSITIONS else "Tv"
            changes.append(f"{orig}>{mut}({sub_type})")
 #Pairs bases of both sequences using the zip function. If they differ among them, it checks if it is a transition or transversion and classifies it accordingly.            
    if not changes:
        return "Synonymous_WT"
 #Classifies mutation sequence as synonymous if both match content.  
    return ",".join(changes)
 #Joints the identified changes into a string and returns it.    
   

def calculate_ts_tv_probability(nuc_change_str, target_length=168, context_length=2686, mode="plasmid_epPCR"):
#Function that calculates the statistical likelihood of a specific mutation occurring base on a binomial probability model. Target sequence length and total length are hard coded. 
    if mode == "genomic_wildtype":
        p_base = 1.0e-9
        ts_fraction = 0.70
        p_indel_ratio = 0.1
    else:  # plasmid_epPCR
        p_base = 1.0e-3
        ts_fraction = 0.80
        p_indel_ratio = 0.01
#Values for type of mutation selected by user for later use in calculations.        

    tv_fraction = 1.0 - ts_fraction
#Calculatin for the fraction of transversions based on value of transitions. 

    if "Indel" in nuc_change_str:
#If the mutation has been classified as a indel, it calculates the probabilit of exactly one insertion ocurring across the target_length.    
   
        prob = math.comb(target_length, 1) * ((p_base * p_indel_ratio) ** 1) * ((1.0 - p_base) ** (target_length - 1))
    elif "Synonymous_WT" in nuc_change_str:
        prob = (1.0 - p_base) ** target_length
#If it is a synonymous sequence, itt calculates the probability of zero changes across the length of the sequence.    
   
    else:
        changes = nuc_change_str.split(",")       
        num_ts = sum(1 for c in changes if "(Ts)" in c)
        num_tv = sum(1 for c in changes if "(Tv)" in c)
        k_total = num_ts + num_tv
#With the split function, it counts and stores in variables the number of transitions and transversions for further calculations.          

        p_ts = p_base * ts_fraction
        p_tv = p_base * tv_fraction
#It takes the p_base value depending on the mode selected by user and multiplying by the fraction of transitions and transversion to get a value of a probability of any base to mutate into a transition or transversion.       

        n_comb = math.comb(target_length, k_total)
        prob = n_comb * (p_ts ** num_ts) * (p_tv ** num_tv) * ((1.0 - p_base) ** (target_length - k_total))
#This applies the binomial probability distribution by calculating n_comb(different ways the total mutaionts could be distributed physcally accross the entire target length), the probability of the transitions and transversions occurring,  
#and finally the probability of other mutatiosn occurring in the rest of the DNA sequence.       

    target_fraction = target_length / context_length    
    return prob * target_fraction
#Calculates the probability obtained by the fraction of the length of sequence allowed to mutate with the total genetic sequence. And returns the final value of the calculation