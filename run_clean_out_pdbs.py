#!/usr/bin/python
__author__ = "morganlnance"


###########################
#### PROGRAM ARGUMENTS ####
###########################

import argparse
parser = argparse.ArgumentParser(description="Use Python to count contacts.")
parser.add_argument("pdb_name_list", help="a file of the pdbs to be analyzed")
parser.add_argument("--ignore_glycosylated_proteins", "-i", action="store_true", help="do you want to skip PDBs that have a covalently attached HETATM group? This is most likely a glycan")
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
import clean_out_pdbs



######################
#### MAIN PROGRAM ####
######################

def go( pdb_name_list, ignore_glycosylated_proteins, heavy_atoms, download_pdbs, keep_cifs, keep_pdbs, keep_clean_pdbs ):
    # get any filenames in the 'pdbs' directory as to not delete them later
    keep_these_files = os.listdir( "pdbs" )
    
    # relevant variable instatiations
    all_clean_pdb_names = []
    
    # load up PDB name list from input argument
    pdb_names_from_list = []
    pdb_file_CR = open( pdb_name_list, 'r' ).readlines()  # CR = carriage return (new line character)
    for name in pdb_file_CR:
        name = name.rstrip( '\n' )
        # add .pdb extension if not already there  -  don't need to do this if downloading the PDB
        if not download_pdbs:
            if not name.endswith( ".pdb" ):
                name = name + ".pdb"
        pdb_names_from_list.append( name )
    
    # get current files as to not delete them later
    working_dir = os.getcwd() + '/'
    cur_files_working = os.listdir( working_dir )
    
    # holds the name of PDBs that couldn't be downloaded
    unable_to_download_pdb_names = []
    
    # instantiate an instance of the Clean object
    # this also instantiates the data holders for PDBs that don't pass the qualifications
    clean = clean_out_pdbs.Clean()
    
    # for each PDB in the list, run the contact counter
    for pdb in pdb_names_from_list:
        if pdb != '':
            # try to print in color if user has colorama library
            try:
                text = "Working on %s" %pdb
                print(Fore.RED + text + Style.RESET_ALL)
            except:
                print "Working on", pdb
            
            # store pdb name in the contact counting class
            pdb_name = pdb[ 0:4 ]
            
            # check to see if this PDB has already been looked at in this round
            if not pdb_name in all_clean_pdb_names:
                # also check lower case
                if not pdb_name.lower() in all_clean_pdb_names:

                    # download the pdb if needed
                    if download_pdbs:
                        pdb = clean.download_pdb( pdb )
                            
                    # check to see if the PDB path exists, otherwise it needed to be downlaoded
                    if not os.path.isfile( pdb ):
                        print "## Skipping", pdb.split( '/' )[-1][0:4], "because it doesn't seem to exist"
                        unable_to_download_pdb_names.append( pdb_name )
                        pass
                        
                    # otherwise the PDB exists and the program continues
                    else:
                        # instantiate data holders for PDB information
                        clean.instantiate_pdb_info_holders()
                        
                        # split ATOM, HETATM, and LINK lines
                        response = clean.split_pdb_file( pdb, ignore_glycosylated_proteins )
                        
                        # if splitting the pdb was successful ( there is a ligand, no AA as ligand, no metal as ligand, no UNK residues, etc. )
                        if response:
                            # get ligand residue numbers from pose
                            response = clean.get_ligand_residues( heavy_atoms, pdb )
                            
                            # if it got this far, then the PDB has an acceptable ligand
                            if response:
                                # write out a clean PDB file
                                if keep_clean_pdbs:
                                    clean.write_clean_pdb_file( pdb )
                                
                                # add the PDB name to the list of clean PDBs
                                all_clean_pdb_names.append( pdb.split( '/' )[-1][0:4] )
                                    
        # get file names in the 'pdbs' directory
        os.chdir( working_dir + 'pdbs' )
        pdb_files = os.listdir( os.getcwd() )
        
        # delete undownloaded .tar.gz files
        for f in pdb_files:
            if f.endswith( ".gz" ):
                os.remove( f )

        # delete unwanted .cif, .pdb, and .pdb files
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
        
        # return to working directory
        os.chdir( working_dir )

    # write out all_clean_pdb_names to a file as those were the PDBs actually analyzed
    with open( "clean_PSMDB_90_non_red_pro_70_non_red_lig_13_ha_cutoff_list", 'wb' ) as fh:
        for line in all_clean_pdb_names:
            fh.write( line + '\n' )
            

 
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
        for line in clean.glycosylated_proteins:
            fh.write( line + '\n' )
        fh.write( "\n" )
        
        fh.write( "## PDBs with an unknown residues\n" )
        for line in clean.unknown_res_pdb_names:
            fh.write( line + '\n' )
        fh.write( "\n" )

        fh.write( "## PDBs with deuterium\n" )
        for line in clean.deuterium_pdb_names:
            fh.write( line + '\n' )
        fh.write( "\n" )

        fh.write( "## PDBs with multiple MODELs\n" )
        for line in clean.multiple_models_pdb_names:
            fh.write( line + '\n' )
        
    # write out the unique three letter codes for ligands that were kept
    filename = "ligand_residues_that_were_kept.txt"
    with open( filename, 'wb' ) as fh:
        for line in clean.uniq_lig_three_letter_codes_kept:
            fh.write( line + '\n' )
        
    # write out the unique three letter codes for ligands that were skipped due to size or linkage
    filename = "ligand_residues_that_were_skipped.txt"
    with open( filename, 'wb' ) as fh:
        fh.write( "## ligand residues that were smaller than %s heavy atoms\n" %heavy_atoms )
        for line in clean.uniq_lig_three_letter_codes_skipped_size:
            fh.write( line + '\n' )
        fh.write( '\n' )
        
        fh.write( "## ligand residues that were skipped because they were linked to the protein\n" )
        for line in clean.uniq_lig_three_letter_codes_skipped_link:
            fh.write( line + '\n' )
    
    print "\nDone! :D\n"




######################
#### RUNS PROGRAM ####
######################

go( input_args.pdb_name_list, input_args.ignore_glycosylated_proteins, input_args.heavy_atoms, input_args.download_pdbs, input_args.keep_cifs, input_args.keep_pdbs, input_args.keep_clean_pdbs )
