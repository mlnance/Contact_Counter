#!/usr/bin/python
__author__ = "morganlnance"


###########################
#### PROGRAM ARGUMENTS ####
###########################

import argparse
parser = argparse.ArgumentParser(description="Use Python to analyze a protein's activesite.")
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
import non_rosetta_analyze_activesite as AS_code



######################
#### MAIN PROGRAM ####
######################

def go( pdb_name_list, ignore_glycosylated_proteins, cutoff, heavy_atoms, download_pdbs, keep_cifs, keep_pdbs, keep_clean_pdbs ):
    # get file names that are already in the 'pdbs' directory (as to not delete them later)
    keep_these_files = os.listdir( "pdbs" )
    
    # relevant variable instatiations
    all_pdb_names = []
    AS = AS_code.ACTIVESITE( pdb_name_list, download_pdbs )
    
    # get current files as to not delete them later
    working_dir = os.getcwd() + '/'
    cur_files_working = os.listdir( working_dir )
    
    # holds the name of PDBs that couldn't be downloaded
    unable_to_download_pdb_names = []
    
    # instantiate the data holders that will hold the data for all of the PDBs
    AS.instantiate_data_holders()
    
    # for each PDB in the list, run the contact counter
    for pdb in AS.pdb_names:
        if pdb != '':
            # try to print in color if user has colorama library
            try:
                text = "Working on %s" %pdb
                print(Fore.RED + text + Style.RESET_ALL)
            except:
                print "Working on", pdb
            
            # instantiate holders for the HETATM and ATOM lines
            AS.instantiate_holders()
            
            # store pdb name in the contact counting class
            pdb_name = pdb[ 0:4 ]
            AS.name = pdb_name
            
            # check to see if this PDB has already been looked at in this round
            if not pdb_name in all_pdb_names:
                # also check lower case
                if not pdb_name.lower() in all_pdb_names:

                    # download the pdb if needed
                    if download_pdbs:
                        pdb = AS.download_pdb( pdb )
                            
                    # check to see if the PDB path exists, otherwise it needed to be downlaoded
                    if not os.path.isfile( pdb ):
                        print "## Skipping", pdb.split( '/' )[-1][0:4], "because it doesn't seem to exist"
                        unable_to_download_pdb_names.append( pdb_name )
                        pass
                        
                    # otherwise the PDB exists and the program continues
                    else:
                        # split ATOM, HETATM, and LINK lines
                        response = AS.split_pdb_file( pdb, ignore_glycosylated_proteins )
                        
                        # if splitting the pdb was successful ( there is a ligand, no AA as ligand, no metal as ligand, no UNK residues, etc. )
                        if response:
                            # get ligand residue numbers from pose
                            response = AS.get_ligand_residues( heavy_atoms, cutoff, pdb, keep_clean_pdbs )
                            
                            # if a ligand remains after the heavy atom cutoff
                            if response:
                                # get protein atoms in the activesite around the ligand
                                response = AS.get_activesite( cutoff )
                                
                                # if there is indeed an activesite
                                if response:
                                    # get activesite composition
                                    AS.get_activesite_AA_composition()
                                    
                                    # get activesite composition per ligandresidue
                                    AS.get_activesite_AA_composition_per_lig_res()
                                    
                                    # count contacts
                                    AS.count_contacts( cutoff )
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
        for line in AS.glycosylated_proteins:
            fh.write( line + '\n' )
        fh.write( "\n" )
    
        fh.write( "## PDBs with an unknown residues\n" )
        for line in AS.unknown_res_pdb_names:
            fh.write( line + '\n' )
        fh.write( "\n" )

        fh.write( "## PDBs with deuterium\n" )
        for line in AS.deuterium_pdb_names:
            fh.write( line + '\n' )
        fh.write( "\n" )

        fh.write( "## PDBs with multiple MODELs\n" )
        for line in AS.multiple_models_pdb_names:
            fh.write( line + '\n' )
    
    # filename will be taken from the name of the PDB list passed through
    if pdb_name_list.endswith( "list" ):
        filename = pdb_name_list.split( "list" )[0]
    else:
        filename = "program_input_"

    # collect AS composition data in a pandas dataframe
    AS_df = pd.DataFrame()
    AS_df["PDB"] = AS.AS_pdb_names
    AS_df["num_lig_res"] = AS.AS_lig_res
    AS_df["num_lig_atoms"] = AS.AS_lig_atms
    AS_df["num_lig_nonpolar_atoms"] = AS.AS_num_ligand_nonpolar_atoms 
    AS_df["num_lig_polar_atoms"] = AS.AS_num_ligand_polar_atoms
    AS_df["num_lig_metal_atoms"] = AS.AS_num_ligand_metal_atoms
    AS_df["num_lig_unk_atom_type"] = AS.AS_num_ligand_unk_atom_types
    AS_df["num_activesite_res"] = AS.AS_activesite_res
    AS_df["num_activesite_atoms"] = AS.AS_activesite_atms
    AS_df["num_activesite_nonpolar_atoms"] = AS.AS_num_activesite_nonpolar_atoms
    AS_df["num_activesite_polar_atoms"] = AS.AS_num_activesite_polar_atoms
    AS_df["num_activesite_unk_atom_type"] = AS.AS_num_activesite_unk_atom_types
    AS_df["ALA"] = AS.ALA
    AS_df["ALA_%_activesite_composition"] = AS.percentage_activesite_ALA
    AS_df["ALA_%_in_full_protein"] = AS.percentage_tot_ALA
    AS_df["CYS"] = AS.CYS
    AS_df["CYS_%_activesite_composition"] = AS.percentage_activesite_CYS
    AS_df["CYS_%_in_full_protein"] = AS.percentage_tot_CYS
    AS_df["ASP"] = AS.ASP
    AS_df["ASP_%_activesite_composition"] = AS.percentage_activesite_ASP
    AS_df["ASP_%_in_full_protein"] = AS.percentage_tot_ASP
    AS_df["GLU"] = AS.GLU
    AS_df["GLU_%_activesite_composition"] = AS.percentage_activesite_GLU
    AS_df["GLU_%_in_full_protein"] = AS.percentage_tot_GLU
    AS_df["PHE"] = AS.PHE
    AS_df["PHE_%_activesite_composition"] = AS.percentage_activesite_PHE
    AS_df["PHE_%_in_full_protein"] = AS.percentage_tot_PHE
    AS_df["GLY"] = AS.GLY
    AS_df["GLY_%_activesite_composition"] = AS.percentage_activesite_GLY
    AS_df["GLY_%_in_full_protein"] = AS.percentage_tot_GLY
    AS_df["HIS"] = AS.HIS
    AS_df["HIS_%_activesite_composition"] = AS.percentage_activesite_HIS
    AS_df["HIS_%_in_full_protein"] = AS.percentage_tot_HIS
    AS_df["ILE"] = AS.ILE
    AS_df["ILE_%_activesite_composition"] = AS.percentage_activesite_ILE
    AS_df["ILE_%_in_full_protein"] = AS.percentage_tot_ILE
    AS_df["LYS"] = AS.LYS
    AS_df["LYS_%_activesite_composition"] = AS.percentage_activesite_LYS
    AS_df["LYS_%_in_full_protein"] = AS.percentage_tot_LYS
    AS_df["LEU"] = AS.LEU
    AS_df["LEU_%_activesite_composition"] = AS.percentage_activesite_LEU
    AS_df["LEU_%_in_full_protein"] = AS.percentage_tot_LEU
    AS_df["MET"] = AS.MET
    AS_df["MET_%_activesite_composition"] = AS.percentage_activesite_MET
    AS_df["MET_%_in_full_protein"] = AS.percentage_tot_MET
    AS_df["ASN"] = AS.ASN
    AS_df["ASN_%_activesite_composition"] = AS.percentage_activesite_ASN
    AS_df["ASN_%_in_full_protein"] = AS.percentage_tot_ASN
    AS_df["PRO"] = AS.PRO
    AS_df["PRO_%_activesite_composition"] = AS.percentage_activesite_PRO
    AS_df["PRO_%_in_full_protein"] = AS.percentage_tot_PRO
    AS_df["GLN"] = AS.GLN
    AS_df["GLN_%_activesite_composition"] = AS.percentage_activesite_GLN
    AS_df["GLN_%_in_full_protein"] = AS.percentage_tot_GLN
    AS_df["ARG"] = AS.ARG
    AS_df["ARG_%_activesite_composition"] = AS.percentage_activesite_ARG
    AS_df["ARG_%_in_full_protein"] = AS.percentage_tot_ARG
    AS_df["SER"] = AS.SER
    AS_df["SER_%_activesite_composition"] = AS.percentage_activesite_SER
    AS_df["SER_%_in_full_protein"] = AS.percentage_tot_SER
    AS_df["THR"] = AS.THR
    AS_df["THR_%_activesite_composition"] = AS.percentage_activesite_THR
    AS_df["THR_%_in_full_protein"] = AS.percentage_tot_THR
    AS_df["VAL"] = AS.VAL
    AS_df["VAL_%_activesite_composition"] = AS.percentage_activesite_VAL
    AS_df["VAL_%_in_full_protein"] = AS.percentage_tot_VAL
    AS_df["TRP"] = AS.TRP
    AS_df["TRP_%_activesite_composition"] = AS.percentage_activesite_TRP
    AS_df["TRP_%_in_full_protein"] = AS.percentage_tot_TRP
    AS_df["TYR"] = AS.TYR
    AS_df["TYR_%_activesite_composition"] = AS.percentage_activesite_TYR
    AS_df["TYR_%_in_full_protein"] = AS.percentage_tot_TYR
    
    print AS_df
    print "\n\n\n\n"
    AS_df.to_csv( str( working_dir ) + '/' + filename + "activesite_AA_composition_at_" + str( cutoff ) + "_Ang_cutoff_and_" + str( heavy_atoms ) + "_heavy_atom_ligand.csv", index = 0, index_col = 0 )

    
    
    # collect AS composition per ligand residue data in a pandas dataframe
    AS_per_lig_df = pd.DataFrame()
    AS_per_lig_df["PDB"] = AS.AS_pdb_names_per_lig
    AS_per_lig_df["uniq_lig_res_names"] = AS.AS_lig_uniq_res_names_per_lig
    AS_per_lig_df["lig_res_names"] = AS.AS_lig_res_names_per_lig
    AS_per_lig_df["num_lig_atoms"] = AS.AS_lig_atms_per_lig
    AS_per_lig_df["num_lig_nonpolar_atoms"] = AS.AS_lig_num_ligand_nonpolar_atoms
    AS_per_lig_df["num_lig_polar_atoms"] = AS.AS_lig_num_ligand_polar_atoms
    AS_per_lig_df["num_lig_metal_atoms"] = AS.AS_lig_num_ligand_metal_atoms
    AS_per_lig_df["num_lig_unk_atoms"] = AS.AS_lig_num_ligand_unk_atoms
    AS_per_lig_df["num_activesite_res"] = AS.AS_activesite_res_per_lig
    AS_per_lig_df["num_activesite_atoms"] = AS.AS_activesite_atms_per_lig
    AS_per_lig_df["num_activesite_nonpolar_atoms"] = AS.AS_lig_num_activesite_nonpolar_atoms
    AS_per_lig_df["num_activesite_polar_atoms"] = AS.AS_lig_num_activesite_polar_atoms
    AS_per_lig_df["ALA"] = AS.ALA_per_lig
    AS_per_lig_df["ALA_%_activesite_composition"] = AS.percentage_activesite_per_lig_ALA
    AS_per_lig_df["ALA_%_in_full_protein"] = AS.percentage_tot_per_lig_ALA
    AS_per_lig_df["CYS"] = AS.CYS_per_lig
    AS_per_lig_df["CYS_%_activesite_composition"] = AS.percentage_activesite_per_lig_CYS
    AS_per_lig_df["CYS_%_in_full_protein"] = AS.percentage_tot_per_lig_CYS
    AS_per_lig_df["ASP"] = AS.ASP_per_lig
    AS_per_lig_df["ASP_%_activesite_composition"] = AS.percentage_activesite_per_lig_ASP
    AS_per_lig_df["ASP_%_in_full_protein"] = AS.percentage_tot_per_lig_ASP
    AS_per_lig_df["GLU"] = AS.GLU_per_lig
    AS_per_lig_df["GLU_%_activesite_composition"] = AS.percentage_activesite_per_lig_GLU
    AS_per_lig_df["GLU_%_in_full_protein"] = AS.percentage_tot_per_lig_GLU
    AS_per_lig_df["PHE"] = AS.PHE_per_lig
    AS_per_lig_df["PHE_%_activesite_composition"] = AS.percentage_activesite_per_lig_PHE
    AS_per_lig_df["PHE_%_in_full_protein"] = AS.percentage_tot_per_lig_PHE
    AS_per_lig_df["GLY"] = AS.GLY_per_lig
    AS_per_lig_df["GLY_%_activesite_composition"] = AS.percentage_activesite_per_lig_GLY
    AS_per_lig_df["GLY_%_in_full_protein"] = AS.percentage_tot_per_lig_GLY
    AS_per_lig_df["HIS"] = AS.HIS_per_lig
    AS_per_lig_df["HIS_%_activesite_composition"] = AS.percentage_activesite_per_lig_HIS
    AS_per_lig_df["HIS_%_in_full_protein"] = AS.percentage_tot_per_lig_HIS
    AS_per_lig_df["ILE"] = AS.ILE_per_lig
    AS_per_lig_df["ILE_%_activesite_composition"] = AS.percentage_activesite_per_lig_ILE
    AS_per_lig_df["ILE_%_in_full_protein"] = AS.percentage_tot_per_lig_ILE
    AS_per_lig_df["LYS"] = AS.LYS_per_lig
    AS_per_lig_df["LYS_%_activesite_composition"] = AS.percentage_activesite_per_lig_LYS
    AS_per_lig_df["LYS_%_in_full_protein"] = AS.percentage_tot_per_lig_LYS
    AS_per_lig_df["LEU"] = AS.LEU_per_lig
    AS_per_lig_df["LEU_%_activesite_composition"] = AS.percentage_activesite_per_lig_LEU
    AS_per_lig_df["LEU_%_in_full_protein"] = AS.percentage_tot_per_lig_LEU
    AS_per_lig_df["MET"] = AS.MET_per_lig
    AS_per_lig_df["MET_%_activesite_composition"] = AS.percentage_activesite_per_lig_MET
    AS_per_lig_df["MET_%_in_full_protein"] = AS.percentage_tot_per_lig_MET
    AS_per_lig_df["ASN"] = AS.ASN_per_lig
    AS_per_lig_df["ASN_%_activesite_composition"] = AS.percentage_activesite_per_lig_ASN
    AS_per_lig_df["ASN_%_in_full_protein"] = AS.percentage_tot_per_lig_ASN
    AS_per_lig_df["PRO"] = AS.PRO_per_lig
    AS_per_lig_df["PRO_%_activesite_composition"] = AS.percentage_activesite_per_lig_PRO
    AS_per_lig_df["PRO_%_in_full_protein"] = AS.percentage_tot_per_lig_PRO
    AS_per_lig_df["GLN"] = AS.GLN_per_lig
    AS_per_lig_df["GLN_%_activesite_composition"] = AS.percentage_activesite_per_lig_GLN
    AS_per_lig_df["GLN_%_in_full_protein"] = AS.percentage_tot_per_lig_GLN
    AS_per_lig_df["ARG"] = AS.ARG_per_lig
    AS_per_lig_df["ARG_%_activesite_composition"] = AS.percentage_activesite_per_lig_ARG
    AS_per_lig_df["ARG_%_in_full_protein"] = AS.percentage_tot_per_lig_ARG
    AS_per_lig_df["SER"] = AS.SER_per_lig
    AS_per_lig_df["SER_%_activesite_composition"] = AS.percentage_activesite_per_lig_SER
    AS_per_lig_df["SER_%_in_full_protein"] = AS.percentage_tot_per_lig_SER
    AS_per_lig_df["THR"] = AS.THR_per_lig
    AS_per_lig_df["THR_%_activesite_composition"] = AS.percentage_activesite_per_lig_THR
    AS_per_lig_df["THR_%_in_full_protein"] = AS.percentage_tot_per_lig_THR
    AS_per_lig_df["VAL"] = AS.VAL_per_lig
    AS_per_lig_df["VAL_%_activesite_composition"] = AS.percentage_activesite_per_lig_VAL
    AS_per_lig_df["VAL_%_in_full_protein"] = AS.percentage_tot_per_lig_VAL
    AS_per_lig_df["TRP"] = AS.TRP_per_lig
    AS_per_lig_df["TRP_%_activesite_composition"] = AS.percentage_activesite_per_lig_TRP
    AS_per_lig_df["TRP_%_in_full_protein"] = AS.percentage_tot_per_lig_TRP
    AS_per_lig_df["TYR"] = AS.TYR_per_lig
    AS_per_lig_df["TYR_%_activesite_composition"] = AS.percentage_activesite_per_lig_TYR
    AS_per_lig_df["TYR_%_in_full_protein"] = AS.percentage_tot_per_lig_TYR

    print AS_per_lig_df
    print "\n\n\n\n"
    AS_per_lig_df.to_csv( str( working_dir ) + '/' + filename + "activesite_AA_composition_per_ligand_at_" + str( cutoff ) + "_Ang_cutoff_and_" + str( heavy_atoms ) + "_heavy_atom_ligand.csv", index = 0, index_col = 0 )
            
    
    # contact counting data
    CC_df = pd.DataFrame()
    CC_df["pdb_names"] = AS.CC_pdb_names
    CC_df["num_lig_atoms"] = AS.CC_lig_atms
    CC_df["num_activesite_atms"] = AS.CC_activesite_atms
    CC_df["num_pp_contacts_within_cutoff"] = AS.CC_pp_contacts
    CC_df["num_pn_contacts_within_cutoff"] = AS.CC_pn_contacts
    CC_df["num_np_contacts_within_cutoff"] = AS.CC_np_contacts
    CC_df["num_nn_contacts_within_cutoff"] = AS.CC_nn_contacts
    CC_df["num_unk_contacts"] = AS.CC_unk_contacts
    
    print CC_df
    print "\n\n\n\n"
    CC_df.to_csv( str( working_dir ) + '/' + filename + "contact_counts_" + str( cutoff ) + "_Ang_cutoff_and_" + str( heavy_atoms ) + "_heavy_atom_ligand.csv", index = 0, index_col = 0 )
    
    
    # make data lists to add over course of program for contact counts per lig - will be added to pandas df at end
    CC_per_lig_df = pd.DataFrame()
    CC_per_lig_df["pdb_names"] = AS.CC_per_lig_pdb_names
    CC_per_lig_df["lig_names"] = AS.CC_per_lig_lig_names
    CC_per_lig_df["num_lig_atms"] = AS.CC_per_lig_lig_atms
    CC_per_lig_df["num_activesite_atms"] = AS.CC_per_lig_activesite_atms
    CC_per_lig_df["pp_contacts_within_cutoff"] = AS.CC_per_lig_pp_contacts
    CC_per_lig_df["pn_contacts_within_cutoff"] = AS.CC_per_lig_pn_contacts
    CC_per_lig_df["np_contacts_within_cutoff"] = AS.CC_per_lig_np_contacts
    CC_per_lig_df["nn_contacts_within_cutoff"] = AS.CC_per_lig_nn_contacts
    CC_per_lig_df["num_unk_contacts"] = AS.CC_per_lig_unk_contacts
    
    print CC_per_lig_df
    CC_per_lig_df.to_csv( str( working_dir ) + '/' + filename + "contact_counts_per_lig_res_" + str( cutoff ) + "_Ang_cutoff_and_" + str( heavy_atoms ) + "_heavy_atom_ligand.csv", index = 0, index_col = 0 )




######################
#### RUNS PROGRAM ####
######################

go( input_args.pdb_name_list, input_args.ignore_glycosylated_proteins, input_args.cutoff, input_args.heavy_atoms, input_args.download_pdbs, input_args.keep_cifs, input_args.keep_pdbs, input_args.keep_clean_pdbs )
