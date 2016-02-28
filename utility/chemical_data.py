__author__ = "morganlnance"


#######################
#### CHEMICAL DATA ####
#######################

# list of polar and nonpolar atoms
nonpolar_atoms = [ 'C' ]
polar_atoms = [ 'O', 'N', 'S', 'P', 'F', 'Cl', 'CL', 'Br', 'BR', 'I', 'Se', 'SE', 'B' ]

# list of residues to remove if they are a ligand
# ligand meaning it was defined as HETATM in the PDB
metal_list = [ 'B', 'K', 'V', 'Y', 'W', 'U', "LI", "BE", "NA", "MG", "AL", "SI", "CA", "SC", "TI", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "GE", "AS", "RB", "SR", "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN", "SB", "TE", "CS", "BA", "LU", "HF", "TA", "RE", "OS", "IR", "PT", "AU", "HG", "TI", "PB", "BI", "PO", "AT", "FR", "RA", "LR", "RF", "DB", "SG", "BH", "HS", "MT", "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "AC", "TH", "PA", "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM", "MD", "NO" ]

AA_list = [ "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TRY", "VAL" ]

nucleic_acid_list = [ 'DA', 'A', 'DT', 'T', 'DG', 'G', 'DC', 'C' ]
