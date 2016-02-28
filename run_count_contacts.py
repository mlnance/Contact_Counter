#!/usr/bin/python
__author__ = "morganlnance"


###########################
#### PROGRAM ARGUMENTS ####
###########################

import argparse
parser = argparse.ArgumentParser(description="Use Python to count contacts.")
parser.add_argument("pdb_name_list", help="a file of the pdbs to be analyzed")
parser.add_argument("--ignore_glycosylated_proteins", "-i", action="store_true", help="do you want to skip PDBs that have a covalently attached HETATM group? This is most likely a glycan")
parser.add_argument("--cutoff", "-c", type=int, default=5, help="how big do you want the activesite cutoff to be, in angstroms? default = 5")
parser.add_argument("--heavy_atoms", "-ha", type=int, default=10, help="how many heavy atoms does a HETATM residue need to be considered a ligand? default = 10")
parser.add_argument("--download_pdbs", "-d", action="store_true", help="do you need to download the pdbs from the database?")
parser.add_argument("--keep_cifs", "-kc", action="store_true", help="do you want to keep the cif files you download?")
parser.add_argument("--keep_pdbs", "-kp", action="store_true", help="do you want to keep the pdbs you download?")
parser.add_argument("--keep_clean_pdbs", "-kcp", action="store_true", help="do you want to keep the cleaned-up version of the pdbs you are working with?")
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

def go( pdb_name_list, ignore_glycosylated_proteins, cutoff, heavy_atoms, download_pdbs, keep_cifs, keep_pdbs, keep_clean_pdbs ):
    # get file names that are already in the 'pdbs' directory (as to not delete them later)
    keep_these_files = os.listdir( "pdbs" )
    
    # relevant variable instatiations
    all_pdb_names = []
    ctct = contact.CTCT( pdb_name_list, download_pdbs )
    
    # get current files as to not delete them later
    working_dir = os.getcwd() + '/'
    cur_files_working = os.listdir( working_dir )
    
    # holds the name of PDBs that couldn't be downloaded
    unable_to_download_pdb_names = []
    
    # instantiate the data holders that will hold the data for all of the PDBs
    ctct.instantiate_data_holders()
    
    # for each PDB in the list, run the contact counter
    for pdb in ctct.pdb_names:
        if pdb != '':
            # try to print in color if user has colorama library
            try:
                text = "Working on %s" %pdb
                print(Fore.RED + text + Style.RESET_ALL)
            except:
                print "Working on", pdb
            
            # instantiate holders for the HETATM and ATOM lines
            ctct.instantiate_holders()
            
            # store pdb name in the contact counting class
            pdb_name = pdb[ 0:4 ]
            ctct.name = pdb_name
            
            # check to see if this PDB has already been looked at in this round
            if not pdb_name in all_pdb_names:
                # also check lower case
                if not pdb_name.lower() in all_pdb_names:

                    # download the pdb if needed
                    if download_pdbs:
                        pdb = ctct.download_pdb( pdb )
                            
                    # check to see if the PDB path exists, otherwise it needed to be downlaoded
                    if not os.path.isfile( pdb ):
                        print "## Skipping", pdb.split( '/' )[-1][0:4], "because it doesn't seem to exist"
                        unable_to_download_pdb_names.append( pdb_name )
                        pass
                        
                    # otherwise the PDB exists and the program continues
                    else:
                        # split ATOM, HETATM, and LINK lines
                        response = ctct.split_pdb_file( pdb, ignore_glycosylated_proteins )
                        
                        # if splitting the pdb was successful ( there is a ligand, no AA as ligand, no metal as ligand, no UNK residues, etc. )
                        if response:
                            # get ligand residue numbers from pose
                            response = ctct.get_ligand_residues( heavy_atoms, cutoff, pdb, keep_clean_pdbs )
                            
                            # if a ligand remains after the heavy atom cutoff
                            if response:
                                # get protein atoms in the activesite around the ligand
                                response = ctct.get_activesite( cutoff )
                                
                                # if there is indeed an activesite
                                if response:
                                    # count contacts
                                    ctct.count_contacts( cutoff )
                                    all_pdb_names.append( pdb.split( '/' )[-1][0:4] )
                                    
        # get file names in the 'pdbs' directory
        os.chdir( working_dir + 'pdbs' )
        pdb_files = os.listdir( os.getcwd() )
        
        # delete undownloaded .tar.gz files
        for f in pdb_files:
            if f.endswith( ".gz" ):
                os.remove( f )

        # delete unwanted .cif, .pdb, and .clean.pdb files
        if not keep_cifs:
            for f in pdb_files:
                if f not in keep_these_files:
                    if f.endswith( ".cif" ):
                        os.remove( f )
        if not keep_pdbs:
            for f in pdb_files:
                if f not in keep_these_files:
                    if f.endswith( ".pdb" ) and len( f ) == 8:
                        os.remove( f )
        if not keep_clean_pdbs:
            for f in pdb_files:
                if f not in keep_these_files:
                    if f.endswith( ".clean.pdb" ):
                        os.remove( f )
        
        # return to working directory
        os.chdir( working_dir )

    # write out all_pdb_names to a file as those were the PDBs actually analyzed
    with open( "clean_PSMDB_90_non_red_pro_70_non_red_lig_13_ha_cutoff_list", 'wb' ) as fh:
        for line in all_pdb_names:
            fh.write( line )
            fh.write( "\n" )
            
    print "\n\n\n"


 
###########################
##### DATA COLLECTION #####
###########################
    
    # write all the PDB names that were either skipped or unable to be downloaded because of the flags given to the program
    # example) skip glycosylated proteins, skip nucleic acids as ligands, skip PDBs with multiple models
    filename = "PDB_files_that_were_skipped.txt"
    with open( filename, 'wb' ) as fh:
        fh.write( "## PDBs unable to be downloaded\n" )
        for line in unable_to_download_pdb_names:
            fh.write( line + '\n' )
        fh.write( "\n" )

        fh.write( "## PDBs with covalently attached HETATMs\n" )
        for line in ctct.glycosylated_proteins:
            fh.write( line + '\n' )
        fh.write( "\n" )
    
        fh.write( "## PDBs with an unknown residues\n" )
        for line in ctct.unknown_res_pdb_names:
            fh.write( line + '\n' )
        fh.write( "\n" )

        fh.write( "## PDBs with deuterium\n" )
        for line in ctct.deuterium_pdb_names:
            fh.write( line + '\n' )
        fh.write( "\n" )

        fh.write( "## PDBs with multiple MODELs\n" )
        for line in ctct.multiple_models_pdb_names:
            fh.write( line + '\n' )
    
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

go( input_args.pdb_name_list, input_args.ignore_glycosylated_proteins, input_args.cutoff, input_args.heavy_atoms, input_args.download_pdbs, input_args.keep_cifs, input_args.keep_pdbs, input_args.keep_clean_pdbs )
