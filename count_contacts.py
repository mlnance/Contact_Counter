import argparse

parser = argparse.ArgumentParser(description="Use PyRosetta to count contacts.")
parser.add_argument("pdb_name_list", help="a file of the pdbs to be analyzed")
parser.add_argument("-cutoff", type=int, default=5, help="how big do you want the activesite cutoff to be, in angstroms? default = 8")
parser.add_argument("-heavy_atoms", type=int, default=7, help="how many heavy atoms does a HETATM residue need to be considered a ligand? default = 7")
parser.add_argument("--download_pdbs", "-d", action="store_true", help="do you need to download the pdbs from the database?")
parser.add_argument("--keep_pdbs", action="store_true", help="do you want to keep the pdbs you download?")
input_args = parser.parse_args()

print "Loading dependencies..."

# import Rosetta modules
from rosetta import pose_from_pdb, get_fa_scorefxn, \
    standard_packer_task, change_cys_state, \
    Pose, MoveMap, RotamerTrialsMover, MinMover, \
    PyMOL_Mover, AtomID, aa_from_oneletter_code
from toolbox import get_hbonds



# path to molfile_to_params file
# can't get this to import for some reason :(
mol2par_dir = "/Users/Research/Rosetta/main/source/src/python/apps/public/"


# path to pymol executable
# must change for each machine!
pymol_exe_dir = "/Applications/MacPyMOL.app/Contents/MacOS/MacPyMOL"


# other imports
import sys
import os
import shutil  # for moving files
import pandas as pd


# list of polar and nonpolar atoms
nonpolar_atoms = ['C']
polar_atoms = ['O', 'N', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'Se', 'B']

# list of residues to remove if they are a ligand
metal_list = ['B', 'K', 'V', 'Y', 'W', 'U', "LI", "BE", "NA", "MG", "AL", "SI", "CA", "SC", "TI", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "GE", "AS", "RB", "SR", "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN", "SB", "TE", "CS", "BA", "LU", "HF", "TA", "RE", "OS", "IR", "PT", "AU", "HG", "TI", "PB", "BI", "PO", "AT", "FR", "RA", "LR", "RF", "DB", "SG", "BH", "HS", "MT", "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "AC", "TH", "PA", "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM", "MD", "NO"]

AA_list = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TRY", "VAL"]



# define exception codes
exception_codes = { -1: "Ligand is an amino acid", -2: "There was an unknown residues", -3: "There is no ligand", -4: "There were no protein atoms within the cutoff distance of the ligand" }


class CTCT():
    def __init__(self):
        # store current working directory
        self.working_dir = os.getcwd()
        
        # load up pdb name list
        self.pdb_names = []
        pdb_file_CR = open(input_args.pdb_name_list, 'r').readlines()  # CR = carriage return
        for name in pdb_file_CR:
            name = name.rstrip('\n')
            # add .pdb extension if not already there
            if not input_args.download_pdbs:
                if not name[:-3] == ".pdb":
                    name = name + ".pdb"
            self.pdb_names.append( name )

        # instantiate a list that will hold files to be deleted after their usage
        self.files_to_delete = []
            


    def download_pdb(self, pdb_name):
        from util import download_pdb
        
        download_pdb( pdb_name, os.getcwd() )



    def pymol_clean(self, pdb_name):
        # use pymol to remove water residues and add hydrogens
        pymol_command = "%s -cqd 'load %s; remove resn hoh; h_add %s; save %s'" %( pymol_exe_dir, pdb_name, pdb_name[0:4], pdb_name )
        
        os.popen( pymol_command )
        


    def split_pdb_file(self, pdb_name):
        # instantiate lists that will hold relevant PDB data
        self.protein = []
        self.ligand = []
        self.lig_names = []  # list of LG1, LG2, LG3, etc
        self.uniq_ligand = []  # list of uniq lig res name and residue stored in the following manner "resname_reschain_resnum"
        AA_lig = []
        water = []
        metals = []
        unknown = []
        
        # move all ATOM and HETATM lines to their appropriate place
        with open( pdb_name, 'r' ) as pdb_fh:
            pdb = pdb_fh.readlines()
            
            for line in pdb:
                if line[0:4] == "ATOM":
                    self.protein.append( line )
                if line[0:5] == "TER  ":
                    self.protein.append( line )
                
                if line[0:6] == "HETATM":
                    # get residue name and remove whitespace
                    res_name = line[17:20].replace( ' ', '' )
                    
                    # check what each residue is by its name (water, metal, ligand)
                    if res_name == "HOH":
                        water.append( line )
                    elif res_name in metal_list:
                        metals.append( line )
                    elif res_name in AA_list:
                        AA_lig.append( line )
                    elif res_name == "UNK":
                        unknown.append( line )
                    else:
                        res_chain = line[21:22]
                        res_num = line[22:26]
                        uniq_lig = res_name + "_" + res_chain + "_" + res_num
                        
                        if uniq_lig not in self.uniq_ligand:
                            self.uniq_ligand.append( uniq_lig )
                        
                        # replace the name of the uniq ligand with LG1, LG2, LG3, etc
                        #new_lig_name = "LG" + str( self.uniq_ligand.index( uniq_lig ) + 1 )
                        new_lig_name = "LG" + str( self.uniq_ligand.index( uniq_lig ) + 1 )
                        temp_line = list( line )
                        temp_line[17:20] = new_lig_name
                        line = "".join( temp_line )
                        self.ligand.append( line )
                        
                        # add new lig name to list, if not already there
                        if new_lig_name not in self.lig_names:
                            self.lig_names.append( new_lig_name )

        # if there are any protein ligands, return an exception code
        if len( AA_lig ) != 0:
            return -1
        
        # if there is an unknown residues, return an exception code
        if len( unknown ) != 0:
            return -2

        # if there is no ligand, return an exception code
        if len( self.ligand ) == 0:
            return -3
        
        # dump each unique ligand residue as a seperate file for mol file creation
        for lig_name in self.lig_names:
            lig_file = self.working_dir + '/' + lig_name + ".pdb"
            ligand = []
            
            for line in self.ligand:
                if line[17:20] == lig_name:
                    ligand.append( line )
            with open( lig_file, 'wb' ) as fh:
                fh.writelines( line for line in ligand )
            
        return 1

    

    def lig_pdb_to_mol(self):
        # must already have LG<num>.pdb in self.working_dir
        # uses babel command in os.popen command
        for lig in self.lig_names:
            ligand_file = self.working_dir + '/' + lig + ".pdb"
            ligand_mol = self.working_dir + '/' + lig + ".mol"
            
            babel_cmd = "babel -ipdb %s -omol %s -h" %( ligand_file, ligand_mol )
            os.popen( babel_cmd )
            
            # get rid of LG<num>.pdb after because it is no longer necessary
            os.remove( ligand_file )
            


    def lig_mol_to_params(self):
        # this is also where the new pdb files for each ligand residue will be collected!
        self.new_ligand_pdbs = []
        
        # have to move into the directory of the program because it doesn't import for some reason...
        os.chdir( mol2par_dir )
        cur_files = os.listdir( mol2par_dir )  # save all current files so they won't get deleted after
        
        # make a params file for each unique ligand mol file
        for lig in self.lig_names:
            param_name = lig
            ligand_mol = self.working_dir + '/' + lig + ".mol"
            
            mol2par_cmd = "python molfile_to_params.py %s --long-names -n %s" %( ligand_mol, param_name )
            os.popen( mol2par_cmd )
        
        for f in os.listdir( mol2par_dir ):
            if f.endswith( ".params" ):
                shutil.move( f, self.working_dir )
            elif f.endswith( ".pdb" ):
                self.new_ligand_pdbs.append( f )
                shutil.move( f, self.working_dir )
            elif f not in cur_files:
                os.remove( f )
                
        # return to working directory
        os.chdir( self.working_dir )
        


    def make_clean_pdb_file(self, pdb_name):
        # make a clean file that will be loaded into Rosetta using the ligand PDBs from params creation
        self.clean_pdb = self.working_dir + '/' + pdb_name[0:4] + ".clean.pdb"
        with open( self.clean_pdb, 'wb' ) as fh:
            # for protein
            fh.writelines( line for line in self.protein )
            
            # for ligand(s)
            for lig_pdb in self.new_ligand_pdbs:
                lig_pdb = self.working_dir + '/' + lig_pdb
                with open( lig_pdb, 'rb' ) as lig_fh:
                    lig = lig_fh.readlines()
                    fh.writelines( line for line in lig )



    def load_pdb(self, pdb_name):
        # store list of ligand params files
        lig_params = []
        for f in os.listdir( self.working_dir ):
            if f.endswith( ".params" ):
                lig_params.append( f )
        
        # extra options string creation
        ext_opts = "-mute basic -mute core -ignore_waters True"
        for param in lig_params:
            ext_opts = ext_opts + " -in:file:extra_res_fa %s" %param
        
        print "Initializing Rosetta with the following options:", ext_opts    
        from rosetta import init
        init(extra_options=ext_opts)
        
        self.pose = pose_from_pdb( pdb_name )
    


    def calc_distance(self, vec1, vec2):
        # takes two lists of xyz coordinates of two atoms
        from math import sqrt, pow

        x1 = vec1[0]
        y1 = vec1[1]
        z1 = vec1[2]
        x2 = vec2[0]
        y2 = vec2[1]
        z2 = vec2[2]

        dist = sqrt( pow(x2 - x1, 2) + pow(y2 - y1, 2) + pow(z2 - z1, 2) )

        return dist



    def get_ligand_residues(self):
        # dictionary: key = lig residue number, value: list of non-hydrogen atom index numbers
        self.lig_dict = {}
        
        for res_num in range( 1, self.pose.total_residue() + 1 ):
            res = self.pose.residue( res_num )
            if res.is_ligand():
                # only keep ligand residues that have the specified number of heavy atoms
                if res.nheavyatoms() >= input_args.heavy_atoms:
                    lig_atom_nums = []
                    num_atms = res.natoms()
                    for atm_num in range(1, num_atms + 1 ):
                        if not res.atom_is_hydrogen( atm_num ):
                            if not res.is_virtual( atm_num ):
                                lig_atom_nums.append( atm_num )
                    self.lig_dict[ res_num ] = lig_atom_nums
        
        # count the number of heavy ligand atoms by looping through the dictionary
        self.num_lig_atms = 0
        for res_num in self.lig_dict:
            self.num_lig_atms += len( self.lig_dict[ res_num ] )
        
        print "  ", self.name, 
        print "has", len( self.lig_dict.keys() ), 
        print "ligand residues with more than", input_args.heavy_atoms, "heavy atoms",
        print "and", self.num_lig_atms, "non-hydrogen atoms"
        
        
        
    def get_activesite_atoms(self):
        # find all residues within cutoff + 3 distance to speed up counting
        close_pro_residues = []
        for lig_res_num in self.lig_dict:
            lig_res = self.pose.residue( lig_res_num )
            # the xyz coordinates of the atom used as this residue's center
            lig_center = lig_res.nbr_atom_xyz()
            for pro_res_num in range( 1, self.pose.total_residue() + 1 ):
                # just skip any residue number that is in the ligand residue dictionary (ie. is a lig res)
                if pro_res_num not in self.lig_dict:
                    if self.pose.residue( pro_res_num ).is_protein():
                        pro_center = self.pose.residue( pro_res_num ).nbr_atom_xyz()
                        if lig_center.distance( pro_center ) <= input_args.cutoff + 3:
                            close_pro_residues.append( pro_res_num )
                        
                        
        # dictionary: key = pro residue number, value: list of atom index numbers within cutoff distance
        self.activesite_dict = {}
        self.num_activesite_atms = []
        
        for lig_res_num in self.lig_dict:
            lig_atm_nums = self.lig_dict[ lig_res_num ]
            for lig_atm_num in lig_atm_nums:
                lig_atm_xyz = self.pose.residue( lig_res_num ).atom( lig_atm_num ).xyz()
                
                for pro_res_num in close_pro_residues:
                    pro_res = self.pose.residue( pro_res_num )
                    activesite_atom_nums = []
                    num_atms = pro_res.natoms()
                    
                    for pro_atm_num in range( 1, num_atms + 1 ):
                        if not pro_res.atom_is_hydrogen( pro_atm_num ):
                            if not pro_res.is_virtual( pro_atm_num ):
                                pro_atm_xyz = pro_res.atom( pro_atm_num ).xyz()
                                if self.calc_distance( lig_atm_xyz, pro_atm_xyz ) <= input_args.cutoff:
                                    # only keep and count unique atoms
                                    if pro_atm_num not in activesite_atom_nums:
                                        activesite_atom_nums.append( pro_atm_num ) 
                                        
                                    # each atm's xyz is unique, so turn it into a string and store to keep number of unique atoms
                                    if str( pro_atm_xyz ) not in self.num_activesite_atms:
                                        self.num_activesite_atms.append( str( pro_atm_xyz ) )
                    
                    if len( activesite_atom_nums ) != 0:
                        self.activesite_dict[ pro_res_num ] = activesite_atom_nums
                
        if len( self.activesite_dict.keys() ) == 0:
            return -4
        else:
            print "  ", self.name, "has", len( self.activesite_dict.keys() ), "activesite residues",
            print "and", len(self.num_activesite_atms), "non-hydrogen activesite atoms"
            return 1
        
        
        
    def count_contacts(self):
        # must have already found all residues within the activesite
        lig_res_nums = self.lig_dict.keys()
        pro_res_nums = self.activesite_dict.keys()
        
        # ligand to protein!!
        polar_polar = 0
        polar_nonpolar = 0
        nonpolar_polar = 0
        nonpolar_nonpolar = 0
        
        for lig_res_num in lig_res_nums:
            for lig_atm_num in self.lig_dict[ lig_res_num ]:
                lig_atm = self.pose.residue( lig_res_num ).atom_type( lig_atm_num )
                lig_atm_xyz = self.pose.residue( lig_res_num ).atom( lig_atm_num ).xyz()
                
                if lig_atm.element() in nonpolar_atoms:
                    for pro_res_num in pro_res_nums:
                        for pro_atm_num in self.activesite_dict[ pro_res_num ]:
                            pro_atm = self.pose.residue( pro_res_num ).atom_type( pro_atm_num )
                            pro_atm_xyz = self.pose.residue( pro_res_num ).atom( pro_atm_num ).xyz()
                            
                            # double check distances
                            if self.calc_distance( lig_atm_xyz, pro_atm_xyz ) <= input_args.cutoff:
                                
                                if pro_atm.element() in nonpolar_atoms:
                                    nonpolar_nonpolar += 1
                                elif pro_atm.element() in polar_atoms:
                                    nonpolar_polar += 1
                                else:
                                    print pro_atm
                                    print '*' * 20, pro_atm.element(), "is not in any list - please update"
                
                elif lig_atm.element() in polar_atoms or lig_atm.element() in metal_list:
                    for pro_res_num in pro_res_nums:
                        for pro_atm_num in self.activesite_dict[ pro_res_num ]:
                            pro_atm = self.pose.residue( pro_res_num ).atom_type( pro_atm_num )
                            pro_atm_xyz = self.pose.residue( pro_res_num ).atom( pro_atm_num ).xyz()
                            
                            # double check distances
                            if self.calc_distance( lig_atm_xyz, pro_atm_xyz ) <= input_args.cutoff:
                                
                                if pro_atm.element() in nonpolar_atoms:
                                    polar_nonpolar += 1
                                elif pro_atm.element() in polar_atoms:
                                    polar_polar += 1
                                else:
                                    print '*' * 20, pro_atm.element(), "is not in any list - please update"
        
                else:
                    print '*' * 20, lig_atm.element(), "is not in any list - please update"

        
        # store all data in global list
        # because if it got here, all of the data existed
        self.pdb_names_final.append( self.name[0:4] )
        self.lig_atms_final.append( self.num_lig_atms )
        self.activesite_atms_final.append( len( self.num_activesite_atms ) )
        self.pp_contacts.append( polar_polar )
        self.pn_contacts.append( polar_nonpolar )
        self.np_contacts.append( nonpolar_polar )
        self.nn_contacts.append( nonpolar_nonpolar )
        
        self.pp_ratio.append( float( polar_polar ) / float( self.num_lig_atms ) )
        self.pn_ratio.append( float( polar_nonpolar ) / float( self.num_lig_atms ) )
        self.np_ratio.append( float( nonpolar_polar ) / float( self.num_lig_atms ) )
        self.nn_ratio.append( float( nonpolar_nonpolar ) / float( self.num_lig_atms ) )
        
        
        
        
    def go(self):
        # make data lists to add over course of program - will be added to pandas df at end
        self.pdb_names_final = []
        self.lig_atms_final = []
        self.activesite_atms_final = []
        self.pp_contacts = []
        self.pn_contacts = []
        self.np_contacts = []
        self.nn_contacts = []
        self.pp_ratio = []
        self.pn_ratio = []
        self.np_ratio = []
        self.nn_ratio = []
        
        cur_files_working = os.listdir( self.working_dir )
        
        for pdb in self.pdb_names:
            if pdb != "":
                print "Loading pdb", pdb
                
                # download the pdb if needed
                if input_args.download_pdbs:
                    self.download_pdb( pdb )
                pdb = pdb[0:4] + ".pdb"
                
                # use pymol to remove waters and add hydrogens for mol file creation
                ##self.pymol_clean( pdb )
                
                # split ATOM and HETATM in order to make params file for ligand
                response = self.split_pdb_file( pdb )
                
                # check to see that no exception was hit in the splitting process
                if response in exception_codes:
                    print exception_codes[ response ]
                    return 0
                
                # turn ligand file into a mol file
                self.lig_pdb_to_mol()
                
                # make a params file for the ligand
                self.lig_mol_to_params()

                # make a clean pdb file
                self.make_clean_pdb_file( pdb )
                
                # load clean pdb into Rosetta through PyRosetta
                clean_pdb = pdb[0:4] + ".clean.pdb"
                self.load_pdb( clean_pdb )
                self.name = self.pose.pdb_info().name()
                print "Working on", self.name
                
                # get ligand residue numbers from pose
                self.get_ligand_residues()
                
                # get protein atoms in the activesite around the ligand
                response = self.get_activesite_atoms()
                if response in exception_codes:
                    print exception_codes[ response ]
                    return 0
                
                # count contacts
                self.count_contacts()
                
                # print a bunch of * to split the info on the screen
                print '*' * 70
                print
                
                for f in os.listdir( self.working_dir ):
                    if f not in cur_files_working:
                        os.remove( f )
            


        # collect data in pandas dataframe
        self.df = pd.DataFrame()
        
        self.df["PDB"] = self.pdb_names_final
        self.df["num_lig_atoms"] = self.lig_atms_final
        self.df["num_activesite_atoms"] = self.activesite_atms_final
        
        self.df["P-P contacts"] = self.pp_contacts
        self.df["P-P ratio"] = self.pp_ratio
        
        self.df["P-NP contacts"] = self.pn_contacts
        self.df["P-NP ratio"] = self.pn_ratio

        self.df["NP-P contacts"] = self.np_contacts
        self.df["NP-P ratio"] = self.np_ratio
        
        self.df["NP-NP contacts"] = self.nn_contacts
        self.df["NP-NP ratio"] = self.nn_ratio
        
        print self.df
