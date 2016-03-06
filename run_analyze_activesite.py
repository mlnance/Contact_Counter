#!/usr/bin/python
__author__ = "morganlnance"


###########################
#### PROGRAM ARGUMENTS ####
###########################

import argparse
parser = argparse.ArgumentParser(description="Use Python to analyze a protein's activesite.")
parser.add_argument("pickle_directory", type=str, help="give me the path to the protein and ligand pickle directory.")
parser.add_argument("cutoff", type=int, default=5, help="how big do you want the activesite cutoff to be, in angstroms? default = 5")
parser.add_argument("heavy_atoms", type=int, default=10, help="how many heavy atoms does a HETATM residue need to be considered a ligand? default = 10")
input_args = parser.parse_args()



#################
#### IMPORTS ####
#################

import os
import sys
try:
    from colorama import Fore, Style
except:
    pass
try:
    import pandas as pd
except ImportError:
    print "Trouble with imports - do you have pandas? Exiting"
    sys.exit()
import non_rosetta_analyze_activesite as AS_code



######################
#### MAIN PROGRAM ####
######################

def go( pickle_dir, cutoff, heavy_atoms ):
    # get the current working directory
    working_dir = os.getcwd()
    
    # check the integrity of the passed pickle_dir
    if not os.path.isdir( pickle_dir ):
        print "You passed me a nonexistent pickle directory. Please check your path:", pickle_dir
        print "Exiting."
        sys.exit()
        
    # add a tailing '/' to the pickle_dir path if needed
    if not pickle_dir.endswith( '/' ):
        pickle_dir += '/'
    
    # collect the PDB names from the protein and ligand pickle files
    os.chdir( pickle_dir )
    all_pdb_names = []
    for pickle_file in os.listdir( os.getcwd() ):
        pdb_name = pickle_file[:4]
        if pdb_name not in all_pdb_names:
            all_pdb_names.append( pdb_name )
    os.chdir( working_dir )
            
    # instantiate an instance of the contact counting code to instantiate class-specific data holders
    activesite = AS_code.ACTIVESITE()
            
    # for each PDB name collected, run the contact counting code using the pro and lig pickles
    for pdb in all_pdb_names:
        # inform the user of the current PDB being worked on
        try:
            # try to print in color if user has colorama library
            text = "Working on %s" %pdb
            print(Fore.RED + text + Style.RESET_ALL)
        except:
            print "Working on", pdb
            
        # store the pickle file paths given the passed information
        pro_pickle_path = pickle_dir + pdb + "_pro.p"
        lig_pickle_path = pickle_dir + pdb + "_lig.p"
        
        # read in the protein and ligand files for this PDB
        activesite.read_pro_lig_pickles( pro_pickle_path, lig_pickle_path )
        
        # look at the ligand and the protein as a whole
        activesite.analyze_protein()
        activesite.analyze_ligand()

        # get protein atoms in the activesite around the ligand
        activesite.get_activesite( cutoff )

        # analyze the activesite
        activesite.get_activesite_atom_composition()
        activesite.get_activesite_atom_composition_per_lig()
        
        #activesite.get_activesite_AA_composition()
        #activesite.get_activesite_AA_composition_per_lig_res()


    print "\n\n\n"


 
###########################
##### DATA COLLECTION #####
###########################
    
    # prefix filename for activesite data
    filename = "program_input_"

    # collect AS composition data in a pandas dataframe
    AS_df = pd.DataFrame()
    AS_df["PDB"] = activesite.AS_pdb_names_list
    AS_df["num_lig_res"] = activesite.AS_lig_res_list
    AS_df["num_lig_atoms"] = activesite.AS_lig_atoms_list
    AS_df["num_lig_nonpolar_atoms"] = activesite.AS_num_ligand_nonpolar_atoms_list
    AS_df["num_lig_polar_atoms"] = activesite.AS_num_ligand_polar_atoms_list
    AS_df["num_lig_metal_atoms"] = activesite.AS_num_ligand_metal_atoms_list
    AS_df["num_lig_unk_atom_type"] = activesite.AS_num_ligand_unk_atom_types_list
    AS_df["num_activesite_res"] = activesite.AS_activesite_res_list
    AS_df["num_activesite_atoms"] = activesite.AS_activesite_atoms_list
    AS_df["num_activesite_nonpolar_atoms"] = activesite.AS_num_activesite_nonpolar_atoms_list
    AS_df["num_activesite_polar_atoms"] = activesite.AS_num_activesite_polar_atoms_list
    AS_df["num_activesite_unk_atom_type"] = activesite.AS_num_activesite_unk_atom_types_list
    AS_df["ALA"] = activesite.ALA_list
    AS_df["ALA_%_activesite_composition"] = activesite.percentage_activesite_ALA_list
    AS_df["ALA_%_in_full_protein"] = activesite.percentage_tot_ALA_list
    AS_df["CYS"] = activesite.CYS_list
    AS_df["CYS_%_activesite_composition"] = activesite.percentage_activesite_CYS_list
    AS_df["CYS_%_in_full_protein"] = activesite.percentage_tot_CYS_list
    AS_df["ASP"] = activesite.ASP_list
    AS_df["ASP_%_activesite_composition"] = activesite.percentage_activesite_ASP_list
    AS_df["ASP_%_in_full_protein"] = activesite.percentage_tot_ASP_list
    AS_df["GLU"] = activesite.GLU_list
    AS_df["GLU_%_activesite_composition"] = activesite.percentage_activesite_GLU_list
    AS_df["GLU_%_in_full_protein"] = activesite.percentage_tot_GLU_list
    AS_df["PHE"] = activesite.PHE_list
    AS_df["PHE_%_activesite_composition"] = activesite.percentage_activesite_PHE_list
    AS_df["PHE_%_in_full_protein"] = activesite.percentage_tot_PHE_list
    AS_df["GLY"] = activesite.GLY_list
    AS_df["GLY_%_activesite_composition"] = activesite.percentage_activesite_GLY_list
    AS_df["GLY_%_in_full_protein"] = activesite.percentage_tot_GLY_list
    AS_df["HIS"] = activesite.HIS_list
    AS_df["HIS_%_activesite_composition"] = activesite.percentage_activesite_HIS_list
    AS_df["HIS_%_in_full_protein"] = activesite.percentage_tot_HIS_list
    AS_df["ILE"] = activesite.ILE_list
    AS_df["ILE_%_activesite_composition"] = activesite.percentage_activesite_ILE_list
    AS_df["ILE_%_in_full_protein"] = activesite.percentage_tot_ILE_list
    AS_df["LYS"] = activesite.LYS_list
    AS_df["LYS_%_activesite_composition"] = activesite.percentage_activesite_LYS_list
    AS_df["LYS_%_in_full_protein"] = activesite.percentage_tot_LYS_list
    AS_df["LEU"] = activesite.LEU_list
    AS_df["LEU_%_activesite_composition"] = activesite.percentage_activesite_LEU_list
    AS_df["LEU_%_in_full_protein"] = activesite.percentage_tot_LEU_list
    AS_df["MET"] = activesite.MET_list
    AS_df["MET_%_activesite_composition"] = activesite.percentage_activesite_MET_list
    AS_df["MET_%_in_full_protein"] = activesite.percentage_tot_MET_list
    AS_df["ASN"] = activesite.ASN_list
    AS_df["ASN_%_activesite_composition"] = activesite.percentage_activesite_ASN_list
    AS_df["ASN_%_in_full_protein"] = activesite.percentage_tot_ASN_list
    AS_df["PRO"] = activesite.PRO_list
    AS_df["PRO_%_activesite_composition"] = activesite.percentage_activesite_PRO_list
    AS_df["PRO_%_in_full_protein"] = activesite.percentage_tot_PRO_list
    AS_df["GLN"] = activesite.GLN_list
    AS_df["GLN_%_activesite_composition"] = activesite.percentage_activesite_GLN_list
    AS_df["GLN_%_in_full_protein"] = activesite.percentage_tot_GLN_list
    AS_df["ARG"] = activesite.ARG_list
    AS_df["ARG_%_activesite_composition"] = activesite.percentage_activesite_ARG_list
    AS_df["ARG_%_in_full_protein"] = activesite.percentage_tot_ARG_list
    AS_df["SER"] = activesite.SER_list
    AS_df["SER_%_activesite_composition"] = activesite.percentage_activesite_SER_list
    AS_df["SER_%_in_full_protein"] = activesite.percentage_tot_SER_list
    AS_df["THR"] = activesite.THR_list
    AS_df["THR_%_activesite_composition"] = activesite.percentage_activesite_THR_list
    AS_df["THR_%_in_full_protein"] = activesite.percentage_tot_THR_list
    AS_df["VAL"] = activesite.VAL_list
    AS_df["VAL_%_activesite_composition"] = activesite.percentage_activesite_VAL_list
    AS_df["VAL_%_in_full_protein"] = activesite.percentage_tot_VAL_list
    AS_df["TRP"] = activesite.TRP_list
    AS_df["TRP_%_activesite_composition"] = activesite.percentage_activesite_TRP_list
    AS_df["TRP_%_in_full_protein"] = activesite.percentage_tot_TRP_list
    AS_df["TYR"] = activesite.TYR_list
    AS_df["TYR_%_activesite_composition"] = activesite.percentage_activesite_TYR_list
    AS_df["TYR_%_in_full_protein"] = activesite.percentage_tot_TYR_list
    
    print AS_df
    print "\n\n\n\n"
    AS_df.to_csv( str( working_dir ) + '/' + filename + "activesite_AA_composition_at_" + str( cutoff ) + "_Ang_cutoff_and_" + str( heavy_atoms ) + "_heavy_atom_ligand.csv", index = 0, index_col = 0 )

    
    
    # collect AS composition per ligand residue data in a pandas dataframe
    AS_per_lig_df = pd.DataFrame()
    AS_per_lig_df["PDB"] = activesite.AS_pdb_names_per_lig_list
    AS_per_lig_df["uniq_lig_res_names"] = activesite.AS_lig_uniq_res_names_per_lig_list
    AS_per_lig_df["lig_res_names"] = activesite.AS_lig_res_names_per_lig_list
    AS_per_lig_df["num_lig_atoms"] = activesite.AS_num_lig_atoms_per_lig_list
    AS_per_lig_df["num_lig_nonpolar_atoms"] = activesite.AS_lig_num_ligand_nonpolar_atoms_list
    AS_per_lig_df["num_lig_polar_atoms"] = activesite.AS_lig_num_ligand_polar_atoms_list
    AS_per_lig_df["num_lig_metal_atoms"] = activesite.AS_lig_num_ligand_metal_atoms_list
    AS_per_lig_df["num_lig_unk_atoms"] = activesite.AS_lig_num_ligand_unk_atoms_list
    AS_per_lig_df["num_activesite_res"] = activesite.AS_activesite_res_per_lig_list
    AS_per_lig_df["num_activesite_atoms"] = activesite.AS_activesite_atoms_per_lig_list
    AS_per_lig_df["num_activesite_nonpolar_atoms"] = activesite.AS_lig_num_activesite_nonpolar_atoms_list
    AS_per_lig_df["num_activesite_polar_atoms"] = activesite.AS_lig_num_activesite_polar_atoms_list
    AS_per_lig_df["ALA"] = activesite.ALA_per_lig_list
    AS_per_lig_df["ALA_%_activesite_composition"] = activesite.percentage_activesite_per_lig_ALA_list
    AS_per_lig_df["ALA_%_in_full_protein"] = activesite.percentage_tot_per_lig_ALA_list
    AS_per_lig_df["CYS"] = activesite.CYS_per_lig_list
    AS_per_lig_df["CYS_%_activesite_composition"] = activesite.percentage_activesite_per_lig_CYS_list
    AS_per_lig_df["CYS_%_in_full_protein"] = activesite.percentage_tot_per_lig_CYS_list
    AS_per_lig_df["ASP"] = activesite.ASP_per_lig_list
    AS_per_lig_df["ASP_%_activesite_composition"] = activesite.percentage_activesite_per_lig_ASP_list
    AS_per_lig_df["ASP_%_in_full_protein"] = activesite.percentage_tot_per_lig_ASP_list
    AS_per_lig_df["GLU"] = activesite.GLU_per_lig_list
    AS_per_lig_df["GLU_%_activesite_composition"] = activesite.percentage_activesite_per_lig_GLU_list
    AS_per_lig_df["GLU_%_in_full_protein"] = activesite.percentage_tot_per_lig_GLU_list
    AS_per_lig_df["PHE"] = activesite.PHE_per_lig_list
    AS_per_lig_df["PHE_%_activesite_composition"] = activesite.percentage_activesite_per_lig_PHE_list
    AS_per_lig_df["PHE_%_in_full_protein"] = activesite.percentage_tot_per_lig_PHE_list
    AS_per_lig_df["GLY"] = activesite.GLY_per_lig_list
    AS_per_lig_df["GLY_%_activesite_composition"] = activesite.percentage_activesite_per_lig_GLY_list
    AS_per_lig_df["GLY_%_in_full_protein"] = activesite.percentage_tot_per_lig_GLY_list
    AS_per_lig_df["HIS"] = activesite.HIS_per_lig_list
    AS_per_lig_df["HIS_%_activesite_composition"] = activesite.percentage_activesite_per_lig_HIS_list
    AS_per_lig_df["HIS_%_in_full_protein"] = activesite.percentage_tot_per_lig_HIS_list
    AS_per_lig_df["ILE"] = activesite.ILE_per_lig_list
    AS_per_lig_df["ILE_%_activesite_composition"] = activesite.percentage_activesite_per_lig_ILE_list
    AS_per_lig_df["ILE_%_in_full_protein"] = activesite.percentage_tot_per_lig_ILE_list
    AS_per_lig_df["LYS"] = activesite.LYS_per_lig_list
    AS_per_lig_df["LYS_%_activesite_composition"] = activesite.percentage_activesite_per_lig_LYS_list
    AS_per_lig_df["LYS_%_in_full_protein"] = activesite.percentage_tot_per_lig_LYS_list
    AS_per_lig_df["LEU"] = activesite.LEU_per_lig_list
    AS_per_lig_df["LEU_%_activesite_composition"] = activesite.percentage_activesite_per_lig_LEU_list
    AS_per_lig_df["LEU_%_in_full_protein"] = activesite.percentage_tot_per_lig_LEU_list
    AS_per_lig_df["MET"] = activesite.MET_per_lig_list
    AS_per_lig_df["MET_%_activesite_composition"] = activesite.percentage_activesite_per_lig_MET_list
    AS_per_lig_df["MET_%_in_full_protein"] = activesite.percentage_tot_per_lig_MET_list
    AS_per_lig_df["ASN"] = activesite.ASN_per_lig_list
    AS_per_lig_df["ASN_%_activesite_composition"] = activesite.percentage_activesite_per_lig_ASN_list
    AS_per_lig_df["ASN_%_in_full_protein"] = activesite.percentage_tot_per_lig_ASN_list
    AS_per_lig_df["PRO"] = activesite.PRO_per_lig_list
    AS_per_lig_df["PRO_%_activesite_composition"] = activesite.percentage_activesite_per_lig_PRO_list
    AS_per_lig_df["PRO_%_in_full_protein"] = activesite.percentage_tot_per_lig_PRO_list
    AS_per_lig_df["GLN"] = activesite.GLN_per_lig_list
    AS_per_lig_df["GLN_%_activesite_composition"] = activesite.percentage_activesite_per_lig_GLN_list
    AS_per_lig_df["GLN_%_in_full_protein"] = activesite.percentage_tot_per_lig_GLN_list
    AS_per_lig_df["ARG"] = activesite.ARG_per_lig_list
    AS_per_lig_df["ARG_%_activesite_composition"] = activesite.percentage_activesite_per_lig_ARG_list
    AS_per_lig_df["ARG_%_in_full_protein"] = activesite.percentage_tot_per_lig_ARG_list
    AS_per_lig_df["SER"] = activesite.SER_per_lig_list
    AS_per_lig_df["SER_%_activesite_composition"] = activesite.percentage_activesite_per_lig_SER_list
    AS_per_lig_df["SER_%_in_full_protein"] = activesite.percentage_tot_per_lig_SER_list
    AS_per_lig_df["THR"] = activesite.THR_per_lig_list
    AS_per_lig_df["THR_%_activesite_composition"] = activesite.percentage_activesite_per_lig_THR_list
    AS_per_lig_df["THR_%_in_full_protein"] = activesite.percentage_tot_per_lig_THR_list
    AS_per_lig_df["VAL"] = activesite.VAL_per_lig_list
    AS_per_lig_df["VAL_%_activesite_composition"] = activesite.percentage_activesite_per_lig_VAL_list
    AS_per_lig_df["VAL_%_in_full_protein"] = activesite.percentage_tot_per_lig_VAL_list
    AS_per_lig_df["TRP"] = activesite.TRP_per_lig_list
    AS_per_lig_df["TRP_%_activesite_composition"] = activesite.percentage_activesite_per_lig_TRP_list
    AS_per_lig_df["TRP_%_in_full_protein"] = activesite.percentage_tot_per_lig_TRP_list
    AS_per_lig_df["TYR"] = activesite.TYR_per_lig_list
    AS_per_lig_df["TYR_%_activesite_composition"] = activesite.percentage_activesite_per_lig_TYR_list
    AS_per_lig_df["TYR_%_in_full_protein"] = activesite.percentage_tot_per_lig_TYR_list

    print AS_per_lig_df
    AS_per_lig_df.to_csv( str( working_dir ) + '/' + filename + "activesite_AA_composition_per_ligand_at_" + str( cutoff ) + "_Ang_cutoff_and_" + str( heavy_atoms ) + "_heavy_atom_ligand.csv", index = 0, index_col = 0 )



######################
#### RUNS PROGRAM ####
######################

go( input_args.pickle_directory, input_args.cutoff, input_args.heavy_atoms )
