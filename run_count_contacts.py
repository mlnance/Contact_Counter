#!/usr/bin/python

import os
import sys
try:
    from colorama import Fore, Style
except:
    pass
import non_rosetta_count_contacts as contact
import argparse

parser = argparse.ArgumentParser(description="Use Python to count contacts.")
parser.add_argument("pdb_name_list", help="a file of the pdbs to be analyzed")
parser.add_argument("--ignore_glycosylated_proteins", "-i", action="store_true", help="do you want to skip PDBs that have a covalently attached HETATM group? This is most likely a glycan")
parser.add_argument("--cutoff", "-c", type=int, default=5, help="how big do you want the activesite cutoff to be, in angstroms? default = 5")
parser.add_argument("--heavy_atoms", "-ha", type=int, default=10, help="how many heavy atoms does a HETATM residue need to be considered a ligand? default = 10")
parser.add_argument("--download_pdbs", "-d", action="store_true", help="do you need to download the pdbs from the database?")
parser.add_argument("--keep_pdbs", action="store_true", help="do you want to keep the pdbs you download?")
parser.add_argument("--keep_clean_pdbs", action="store_true", help="do you want to keep the cleaned-up version of the pdbs you are working with?")
input_args = parser.parse_args()



def go( pdb_name_list, ignore_glycosylated_proteins, cutoff, heavy_atoms, download_pdbs, keep_pdbs, keep_clean_pdbs ):
    # relevant variable instatiations
    all_pdb_names = []
    ctct = contact.CTCT( pdb_name_list, download_pdbs )
    
    # get current files as to not delete them later
    cur_dir = os.getcwd() + '/'
    cur_files_working = os.listdir( cur_dir )
    
    pdb_name_list = []  # only the first four letters

    # instantiate the data holders that will hold the data for all of the PDBs
    ctct.instantiate_data_holders()
    
    for pdb in ctct.pdb_names:
        if pdb != '':
            try:
                text = "Working on %s" %pdb
                print(Fore.RED + text + Style.RESET_ALL)
            except:
                print "Working on", pdb
            
            # instantiate holders for the HETATM and ATOM lines
            ctct.instantiate_holders()
            
            # get pdb name
            ctct.name = pdb
            
            # check to see if this PDB has already been looked at in this round
            if not pdb[0:4] in pdb_name_list:
                if not pdb[0:4].lower() in pdb_name_list:  # check lower case
                    pdb_name_list.append( pdb )
                    
                    # download the pdb if needed
                    if download_pdbs:
                        pdb = ctct.download_pdb( pdb )
                            
                    # otherwise, try opening the file in the current directory
                    else:
                        try:
                            with open( pdb, 'r' ) as pdb_fh:
                                pdb_lines = pdb_fh.readlines()
                        except IOError:
                            print pdb, "doesn't exist in this directory. Did you mean to download it? Add the '-d' argument in your command. Exiting."
                            sys.exit()

                    ## use pymol to remove waters and hydrogens
                    #pymol_clean( pdb )
                            
                    # split ATOM and HETATM
                    response = ctct.split_pdb_file( pdb, ignore_glycosylated_proteins, keep_clean_pdbs )
                        
                    # if splitting the pdb was successful ( there is a ligand, no AA as ligand, no metal as ligand, no UNK residues )
                    if response:
                        # get ligand residue numbers from pose
                        response = ctct.get_ligand_residues( heavy_atoms, cutoff )
                            
                    # if a ligand remains after the heavy atom cutoff
                    if response:
                        # get protein atoms in the activesite around the ligand
                        response = ctct.get_activesite( cutoff )
                        
                        # if there is indeed an activesite
                        if response:
                            # get activesite composition
                            ctct.get_activesite_AA_composition()
                            
                            # get activesite composition per ligandresidue
                            ctct.get_activesite_AA_composition_per_lig_res()
                            
                            # count contacts
                            ctct.count_contacts( cutoff )
                            all_pdb_names.append( pdb )
                                
        with open( "clean_PSMDB_90_pro_70_lig_7_atom_cutoff_list", 'wb' ) as fh:
            fh.writelines( all_pdb_names )
            
        # PDB files always end with .pdb and are 8 characters long
        # this program creates XXXX.clean.pdb files, so have to be specific on what gets deleted
        if not keep_pdbs:
            os.chdir( cur_dir + 'pdbs' )
            pdb_files = os.listdir( os.getcwd() )
            for f in pdb_files:
                if len( f ) == 8:
                    os.remove( f )
            os.chdir( cur_dir )


 
###########################
##### DATA COLLECTION #####
###########################
'''
        # print the names of glycosylated proteins to a file including the name of the list passed
        filename = pdb_name_list + "_glycosylated_protein_PDB_names.txt"
        with open( filename, 'wb' ) as fh:
            fh.writelines( glycosylated_proteins )
            fh.write( "\n" )
                
        # filename will be taken from the name of the PDB list passed through
        filename = pdb_name_list.split( "list" )[0]
                                
        # collect AA composition data in a pandas dataframe
        AA_df = pd.DataFrame()
        AA_df["PDB"] = AA_pdb_names
        AA_df["num_lig_res"] = AA_lig_res
        AA_df["num_lig_atoms"] = AA_lig_atms
        AA_df["num_lig_nonpolar_atoms"] = AA_num_ligand_nonpolar_atoms 
        AA_df["num_lig_polar_atoms"] = AA_num_ligand_polar_atoms
        AA_df["num_lig_unk_atom_type"] = AA_num_ligand_unk_atom_types
        AA_df["num_activesite_res"] = AA_activesite_res
        AA_df["num_activesite_atoms"] = AA_activesite_atms
        AA_df["num_activesite_nonpolar_atoms"] = AA_num_activesite_nonpolar_atoms
        AA_df["num_activesite_polar_atoms"] = AA_num_activesite_polar_atoms
        AA_df["num_activesite_unk_atom_type"] = AA_num_activesite_unk_atom_types
        AA_df["ALA"] = ALA
        AA_df["CYS"] = CYS
        AA_df["ASP"] = ASP
        AA_df["GLU"] = GLU
        AA_df["PHE"] = PHE
        AA_df["GLY"] = GLY
        AA_df["HIS"] = HIS
        AA_df["ILE"] = ILE
        AA_df["LYS"] = LYS
        AA_df["LEU"] = LEU
        AA_df["MET"] = MET
        AA_df["ASN"] = ASN
        AA_df["PRO"] = PRO
        AA_df["GLN"] = GLN
        AA_df["ARG"] = ARG
        AA_df["SER"] = SER
        AA_df["THR"] = THR
        AA_df["VAL"] = VAL
        AA_df["TRP"] = TRP
        AA_df["TYR"] = TYR
        
        print AA_df
        AA_df.to_csv( str( working_dir ) + '/' + filename + "activesite_AA_composition_at_" + str( cutoff ) + "_Ang_cutoff_and_" + str( heavy_atoms ) + "_heavy_atom_ligand.csv", index = 0, index_col = 0 )


            
        # collect AA composition per ligand residue data in a pandas dataframe
        AA_per_lig_df = pd.DataFrame()
        AA_per_lig_df["PDB"] = AA_pdb_names_per_lig
        AA_per_lig_df["uniq_lig_res_names"] = AA_lig_uniq_res_names_per_lig
        AA_per_lig_df["lig_res_names"] = AA_lig_res_names_per_lig
        AA_per_lig_df["num_lig_atoms"] = AA_lig_atms_per_lig
        AA_per_lig_df["num_lig_nonpolar_atoms"] = AA_lig_num_ligand_nonpolar_atoms
        AA_per_lig_df["num_lig_polar_atoms"] = AA_lig_num_ligand_polar_atoms
        AA_per_lig_df["num_lig_unk_atoms"] = AA_lig_num_ligand_unk_atoms
        AA_per_lig_df["num_activesite_res"] = AA_activesite_res_per_lig
        AA_per_lig_df["num_activesite_atoms"] = AA_activesite_atms_per_lig
        AA_per_lig_df["num_activesite_nonpolar_atoms"] = AA_lig_num_activesite_nonpolar_atoms
        AA_per_lig_df["num_activesite_polar_atoms"] = AA_lig_num_activesite_polar_atoms
        AA_per_lig_df["ALA"] = ALA_per_lig
        AA_per_lig_df["CYS"] = CYS_per_lig
        AA_per_lig_df["ASP"] = ASP_per_lig
        AA_per_lig_df["GLU"] = GLU_per_lig
        AA_per_lig_df["PHE"] = PHE_per_lig
        AA_per_lig_df["GLY"] = GLY_per_lig
        AA_per_lig_df["HIS"] = HIS_per_lig
        AA_per_lig_df["ILE"] = ILE_per_lig
        AA_per_lig_df["LYS"] = LYS_per_lig
        AA_per_lig_df["LEU"] = LEU_per_lig
        AA_per_lig_df["MET"] = MET_per_lig
        AA_per_lig_df["ASN"] = ASN_per_lig
        AA_per_lig_df["PRO"] = PRO_per_lig
        AA_per_lig_df["GLN"] = GLN_per_lig
        AA_per_lig_df["ARG"] = ARG_per_lig
        AA_per_lig_df["SER"] = SER_per_lig
        AA_per_lig_df["THR"] = THR_per_lig
        AA_per_lig_df["VAL"] = VAL_per_lig
        AA_per_lig_df["TRP"] = TRP_per_lig
        AA_per_lig_df["TYR"] = TYR_per_lig
        
        print AA_per_lig_df
        AA_per_lig_df.to_csv( str( working_dir ) + '/' + filename + "activesite_AA_composition_per_ligand_at_" + str( cutoff ) + "_Ang_cutoff_and_" + str( heavy_atoms ) + "_heavy_atom_ligand.csv", index = 0, index_col = 0 )
            

        # contact counting data
        CC_df = pd.DataFrame()
        CC_df["pdb_names"] = CC_pdb_names
        CC_df["num_lig_atoms"] = CC_lig_atms
        CC_df["num_activesite_atms"] = CC_activesite_atms
        CC_df["num_polar_polar_contacts"] = CC_pp_contacts
        CC_df["num_polar_nonpolar_contacts"] = CC_pn_contacts
        CC_df["num_nonpolar_polar_contacts"] = CC_np_contacts
        CC_df["num_nonpolar_nonpolar_contacts"] = CC_nn_contacts
        CC_df["num_unk_contacts"] = CC_unk_contacts
        
        print CC_df
        CC_df.to_csv( str( working_dir ) + '/' + filename + "contact_counts_" + str( cutoff ) + "_Ang_cutoff_and_" + str( heavy_atoms ) + "_heavy_atom_ligand.csv", index = 0, index_col = 0 )


        # make data lists to add over course of program for contact counts per lig - will be added to pandas df at end
        CC_per_lig_df = pd.DataFrame()
        CC_per_lig_df["pdb_names"] = CC_per_lig_pdb_names
        CC_per_lig_df["lig_names"] = CC_per_lig_lig_names
        CC_per_lig_df["num_lig_atms"] = CC_per_lig_lig_atms
        CC_per_lig_df["num_activesite_atms"] = CC_per_lig_activesite_atms
        CC_per_lig_df["polar_polar_contacts"] = CC_per_lig_pp_contacts
        CC_per_lig_df["polar_nonpolar_contacts"] = CC_per_lig_pn_contacts
        CC_per_lig_df["nonpolar_polar_contacts"] = CC_per_lig_np_contacts
        CC_per_lig_df["nonpolar_nonpolar_contacts"] = CC_per_lig_nn_contacts
        CC_per_lig_df["num_unk_contacts"] = CC_per_lig_unk_contacts

        print CC_per_lig_df
        CC_per_lig_df.to_csv( str( working_dir ) + '/' + filename + "contact_counts_per_lig_res_" + str( cutoff ) + "_Ang_cutoff_and_" + str( heavy_atoms ) + "_heavy_atom_ligand.csv", index = 0, index_col = 0 )
'''



####################
### RUNS PROGRAM ###
####################


#my_obj = CTCT()
go( input_args.pdb_name_list, input_args.ignore_glycosylated_proteins, input_args.cutoff, input_args.heavy_atoms, input_args.download_pdbs, input_args.keep_pdbs, input_args.keep_clean_pdbs )
