"""
Grantham Distance Matrix and Pairwise Scoring
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
    This code calculates the physicochemical differences between normal and mutated protein sequences by applying Grantham´s distance metrics.

================================================================================
"""

GRANTHAM_MATRIX = {
    ('A', 'A'): 0,   ('A', 'R'): 112, ('A', 'N'): 111, ('A', 'D'): 126, ('A', 'C'): 195, 
    ('A', 'E'): 107, ('A', 'Q'): 91,  ('A', 'G'): 60,  ('A', 'H'): 86,  ('A', 'I'): 94,  
    ('A', 'L'): 96,  ('A', 'K'): 106, ('A', 'M'): 84,  ('A', 'F'): 113, ('A', 'P'): 27,  
    ('A', 'S'): 99,  ('A', 'T'): 58,  ('A', 'W'): 148, ('A', 'Y'): 112, ('A', 'V'): 64,
    ('R', 'R'): 0,   ('R', 'N'): 86,  ('R', 'D'): 96,  ('R', 'C'): 180, ('R', 'E'): 54,  
    ('R', 'Q'): 43,  ('R', 'G'): 125, ('R', 'H'): 29,  ('R', 'I'): 97,  ('R', 'L'): 102, 
    ('R', 'K'): 26,  ('R', 'M'): 91,  ('R', 'F'): 97,  ('R', 'P'): 103, ('R', 'S'): 110, 
    ('R', 'T'): 71,  ('R', 'W'): 101, ('R', 'Y'): 77,  ('R', 'V'): 96,
    ('N', 'N'): 0,   ('N', 'D'): 23,  ('N', 'C'): 139, ('N', 'E'): 42,  ('N', 'Q'): 46,  
    ('N', 'G'): 80,  ('N', 'H'): 68,  ('N', 'I'): 149, ('N', 'L'): 153, ('N', 'K'): 94,  
    ('N', 'M'): 142, ('N', 'F'): 158, ('N', 'P'): 91,  ('N', 'S'): 46,  ('N', 'T'): 65,  
    ('N', 'W'): 174, ('N', 'Y'): 143, ('N', 'V'): 133,
    ('D', 'D'): 0,   ('D', 'C'): 154, ('D', 'E'): 45,  ('D', 'Q'): 61,  ('D', 'G'): 108, 
    ('D', 'H'): 81,  ('D', 'I'): 168, ('D', 'L'): 172, ('D', 'K'): 101, ('D', 'M'): 160, 
    ('D', 'F'): 177, ('D', 'P'): 108, ('D', 'S'): 65,  ('D', 'T'): 85,  ('D', 'W'): 181, 
    ('D', 'Y'): 160, ('D', 'V'): 152,
    ('C', 'C'): 0,   ('C', 'E'): 170, ('C', 'Q'): 154, ('C', 'G'): 159, ('C', 'H'): 174, 
    ('C', 'I'): 198, ('C', 'L'): 198, ('C', 'K'): 202, ('C', 'M'): 196, ('C', 'F'): 205, 
    ('C', 'P'): 169, ('C', 'S'): 112, ('C', 'T'): 149, ('C', 'W'): 215, ('C', 'Y'): 194, 
    ('C', 'V'): 192,
    ('E', 'E'): 0,   ('E', 'Q'): 29,  ('E', 'G'): 98,  ('E', 'H'): 40,  ('E', 'I'): 134, 
    ('E', 'L'): 138, ('E', 'K'): 56,  ('E', 'M'): 126, ('E', 'F'): 140, ('E', 'P'): 93,  
    ('E', 'S'): 80,  ('E', 'T'): 65,  ('E', 'W'): 152, ('E', 'Y'): 122, ('E', 'V'): 121,
    ('Q', 'Q'): 0,   ('Q', 'G'): 87,  ('Q', 'H'): 24,  ('Q', 'I'): 109, ('Q', 'L'): 113, 
    ('Q', 'K'): 53,  ('Q', 'M'): 87,  ('Q', 'F'): 116, ('Q', 'P'): 76,  ('Q', 'S'): 68,  
    ('Q', 'T'): 42,  ('Q', 'W'): 130, ('Q', 'Y'): 99,  ('Q', 'V'): 96,
    ('G', 'G'): 0,   ('G', 'H'): 89,  ('G', 'I'): 135, ('G', 'L'): 138, ('G', 'K'): 127, 
    ('G', 'M'): 127, ('G', 'F'): 153, ('G', 'P'): 42,  ('G', 'S'): 56,  ('G', 'T'): 59,  
    ('G', 'W'): 184, ('G', 'Y'): 147, ('G', 'V'): 109,
    ('H', 'H'): 0,   ('H', 'I'): 94,  ('H', 'L'): 99,  ('H', 'K'): 32,  ('H', 'M'): 87,  
    ('H', 'F'): 100, ('H', 'P'): 77,  ('H', 'S'): 89,  ('H', 'T'): 47,  ('H', 'W'): 115, 
    ('H', 'Y'): 83,  ('H', 'V'): 84,
    ('I', 'I'): 0,   ('I', 'L'): 5,   ('I', 'K'): 102, ('I', 'M'): 10,  ('I', 'F'): 21,  
    ('I', 'P'): 95,  ('I', 'S'): 142, ('I', 'T'): 89,  ('I', 'W'): 61,  ('I', 'Y'): 33,  
    ('I', 'V'): 29,
    ('L', 'L'): 0,   ('L', 'K'): 107, ('L', 'M'): 15,  ('L', 'F'): 22,  ('L', 'P'): 98,  
    ('L', 'S'): 145, ('L', 'T'): 92,  ('L', 'W'): 61,  ('L', 'Y'): 36,  ('L', 'V'): 32,
    ('K', 'K'): 0,   ('K', 'M'): 95,  ('K', 'F'): 102, ('K', 'P'): 103, ('K', 'S'): 121, 
    ('K', 'T'): 78,  ('K', 'W'): 110, ('K', 'Y'): 85,  ('K', 'V'): 97,
    ('M', 'M'): 0,   ('M', 'F'): 28,  ('M', 'P'): 87,  ('M', 'S'): 135, ('M', 'T'): 81,  
    ('M', 'W'): 67,  ('M', 'Y'): 36,  ('M', 'V'): 21,
    ('F', 'F'): 0,   ('F', 'P'): 114, ('F', 'S'): 155, ('F', 'T'): 103, ('F', 'W'): 40,  
    ('F', 'Y'): 22,  ('F', 'V'): 50,
    ('P', 'P'): 0,   ('P', 'S'): 74,  ('P', 'T'): 38,  ('P', 'W'): 147, ('P', 'Y'): 115, 
    ('P', 'V'): 68,
    ('S', 'S'): 0,   ('S', 'T'): 58,  ('S', 'W'): 177, ('S', 'Y'): 144, ('S', 'V'): 124,
    ('T', 'T'): 0,   ('T', 'W'): 128, ('T', 'Y'): 92,  ('T', 'V'): 29,
    ('W', 'W'): 0,   ('W', 'Y'): 37,  ('W', 'V'): 88,
    ('Y', 'Y'): 0,   ('Y', 'V'): 55,
    ('V', 'V'): 0
}
#A dictionary that connents the different possible mutations to the Grantham distance value. It only contains changes in one direction. 

def get_grantham_distance(aa1, aa2):
    if aa1 == aa2:
        return 0
    return GRANTHAM_MATRIX.get((aa1, aa2), GRANTHAM_MATRIX.get((aa2, aa1), None))
#This function determines the value of the grantham distance according to the amino acid swap. If inputs are identical, it returns a 0. Otherwise, it gets the result, but because the dictionary only has one direction data, if the specific
#amino acid pair is not found, the function reverse the order to find it. If neither combination exists in the dictionary, it returns a None value. 

def calculate_grantham_score(wt_prot, mut_prot, gap_penalty=100):
    wt_clean = wt_prot.replace("*", "").replace("-", "")
    mut_clean = mut_prot.replace("*", "").replace("-", "")
#First, the function cleaning input by removing asterisks and dashes(stop codons or gaps) that could affect processing. It applies a value of 100 for penalty,as an average score for all possible values. 

    if not mut_clean:
        return len(wt_clean) * gap_penalty
#If sequences ais empty, like a early STOP codon, it multiples the total length of wildtype sequence by the gap penalty. 

    if len(wt_clean) == len(mut_clean):
        total_dist = 0
        for a, b in zip(wt_clean, mut_clean):
            dist = get_grantham_distance(a, b)
            total_dist += dist if dist is not None else gap_penalty
        return total_dist
#When sequences have the same lenght, it compares amino acids through the zip function, calculating Grantham distance for each and adding them up
#to the total_dist. When there is a None character in the sequence, if applyes the gap penalty. 

    min_len = min(len(wt_clean), len(mut_clean))
    total_dist = 0
#If lenghts do not match, it determines the lenght of the shorter sequence in min_len.     
    for i in range(min_len):
        dist = get_grantham_distance(wt_clean[i], mut_clean[i])
        total_dist += dist if dist is not None else gap_penalty
#This section does the iteration and calculation of the Grantham distance in shorter distances to account for earlier STOP codons.    
    
    length_diff = abs(len(wt_clean) - len(mut_clean))
    total_dist += length_diff * gap_penalty
#This section calculates the penalty value taken into account the lenght difference between the wildtype and mutant sequence and adds it to the total_dist 
   
    return total_dist
#Functio returns the calculated Grantham score for the sequence.     