#!/usr/bin/python
__author__ = "morganlnance"


###########################
#### PROGRAM ARGUMENTS ####
###########################

import argparse
parser = argparse.ArgumentParser(description="Use Python to count contacts.")
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
import non_rosetta_count_contacts as contact



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
    ctct = contact.CTCT()
            
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
        ctct.read_pro_lig_pickles( pro_pickle_path, lig_pickle_path )

        # get protein atoms in the activesite around the ligand
        ctct.get_activesite( cutoff )

        # count contacts
        ctct.count_contacts( cutoff )

                                    
    print "\n\n\n"


 
###########################
##### DATA COLLECTION #####
###########################
    
    # filename will be taken from the name of the PDB list passed through
    if pdb_name_list.endswith( "list" ):
        filename = pdb_name_list.split( "list" )[0]
    else:
        filename = "program_input_"

    # contact counting data
    CC_df = pd.DataFrame()
    CC_df["pdb_names"] = ctct.CC_pdb_names
    CC_df["num_lig_atoms"] = ctct.CC_lig_atms
    CC_df["num_activesite_atms"] = ctct.CC_activesite_atms
    CC_df["num_pp_contacts_within_cutoff"] = ctct.CC_pp_contacts
    CC_df["num_pn_contacts_within_cutoff"] = ctct.CC_pn_contacts
    CC_df["num_np_contacts_within_cutoff"] = ctct.CC_np_contacts
    CC_df["num_nn_contacts_within_cutoff"] = ctct.CC_nn_contacts
    CC_df["num_unk_contacts"] = ctct.CC_unk_contacts
    
    print CC_df
    print "\n\n\n\n"
    CC_df.to_csv( str( working_dir ) + '/' + filename + "contact_counts_" + str( cutoff ) + "_Ang_cutoff_and_" + str( heavy_atoms ) + "_heavy_atom_ligand.csv", index = 0, index_col = 0 )
    
    
    # make data lists to add over course of program for contact counts per lig - will be added to pandas df at end
    CC_per_lig_df = pd.DataFrame()
    CC_per_lig_df["pdb_names"] = ctct.CC_per_lig_pdb_names
    CC_per_lig_df["lig_names"] = ctct.CC_per_lig_lig_names
    CC_per_lig_df["num_lig_atms"] = ctct.CC_per_lig_lig_atms
    CC_per_lig_df["num_activesite_atms"] = ctct.CC_per_lig_activesite_atms
    CC_per_lig_df["pp_contacts_within_cutoff"] = ctct.CC_per_lig_pp_contacts
    CC_per_lig_df["pn_contacts_within_cutoff"] = ctct.CC_per_lig_pn_contacts
    CC_per_lig_df["np_contacts_within_cutoff"] = ctct.CC_per_lig_np_contacts
    CC_per_lig_df["nn_contacts_within_cutoff"] = ctct.CC_per_lig_nn_contacts
    CC_per_lig_df["num_unk_contacts"] = ctct.CC_per_lig_unk_contacts
    
    print CC_per_lig_df
    CC_per_lig_df.to_csv( str( working_dir ) + '/' + filename + "contact_counts_per_lig_res_" + str( cutoff ) + "_Ang_cutoff_and_" + str( heavy_atoms ) + "_heavy_atom_ligand.csv", index = 0, index_col = 0 )





######################
#### RUNS PROGRAM ####
######################

go( input_args.pickle_directory, input_args.cutoff, input_args.heavy_atoms )
