# ContactCounter

## Default values for optional arguments
-cutoff:				 activesite cutoff distance  =  5 Angstroms

-heavy_atoms:          			 number of heavy atoms to be considered a ligand  =  10

--download_pdbs, -d:   		       	 download PDBs from internet  =  False

--keep_pdbs:           			 keep PDB files after download  =  False

--keep_clean_pdbs:     		     	 keep cleaned-up versions of the PDB files after use  =  False

--ignore_glycosylated_proteins, -i:	 when determining a ligand, ignore covalenlty bound residues = False


## Example ways to run this program
1) Use default arguments and no downloading of PDBs

./program_name pdb_list

2) Use default arguments and download PDBs

./program_name pdb_list --download_pdbs/-d

3) Don't download PDBs but change activesite cutoff distance to 7 Ang and heavy atom cutoff to 13

./program_name pdb_list -cutoff=7 -heavy_atoms=13

4) Download and keep PDBs using default values

./program_name pdb_list -d --keep_pdbs

5) Download and keep both downloaded and cleaned-up PDBs using default values

./program_name pdb_list -d --keep_pdbs --keep_clean_pdbs


# Function overview
## split_pdb_file
self.protein_lines : an instance of PDB_line for each line in the PDB starting with ATOM
self.ligand_lines : an instance of PDB_line for each line in the PDB starting with HETATM, unless it is a water residue, metal residue, amino acid residue, or an unknown residue
self.protein : a dictionary where each key is a unique protein name (resname_reschain_resnum) and the values are the PDB_line instances of each corresponding line from the PDB
self.ligand : a dictionary where each key is a unique ligand name (resname_reschain_resnum) and the values are the PDB_line instances of each corresponding line from the PDB

## get_ligand_residues
self.ligand_dict : a dictionary where each key is a unqiue ligand name (resname_reschain_resnum) that passes the user specified heavy atom cutoff, and the values are the PDB_line instances of each corresponding line from the PDB
self.lig_res_names : a list of the 3-letter names for each ligand, where repeats are allowed
self.uniq_lig_res_names : a list of the 3-letter names for each ligand, where no repeats are allowed - this is to be used to count the number of each ligand residue type in the protein. Essentially, this is the result of self.lig_res_names.count( "name" )
self.num_ligand_residues : the number of ligand residues that satisfied the heavy atom cutoff in the protein
self.num_ligand_atoms : the number of atoms in each ligand residue from self.num_ligand_residues

## get_activesite
self.activesite_dict : a dictionary where each key is a unique protein name (resname_reschain_resnum) and each value is a list of the corresponding PDB lines
self.activesite_lig_pro_res_dict : a dictionary where each key is a unique ligand name (resname_reschain_resnum) and each value is a list of the 3-letter name of each protein residue within the specified cutoff distance
self.activesite_lig_pro_atms_dict : a dictionary where each key is a unique ligand name (resname_reschain_resnum) and each value is a list of the PDB lines of each protein residue within the specified cutoff distance
self.num_activesite_res : the number of residues within the user cutoff's activesite for each ligand residue
self.num_activesite_atms : the number of atoms of each residue within the user cutoff's activesite for each ligand residue
self.activesite_residues : the unique name (resname_reschain_resnum) of each protein residue within the user cutoff's activesite
