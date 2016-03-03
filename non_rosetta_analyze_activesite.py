#!/usr/bin/python
__author__ = "morganlnance"


#####
# path to pymol executable
# must change for each machine!
#pymol_exe_dir = "/Applications/MacPyMOL.app/Contents/MacOS/MacPyMOL"
#####



#################
#### IMPORTS ####
#################

print "Loading '%s' dependencies...\n" %__name__

import sys
import os
import shutil
import pickle
try:
    from colorama import Fore, Style
except:
    pass

# utility functions
sys.path.append( "utility" )
from line_definitions import *
from chemical_data import *
from util import calc_distance



##################################
#### ANALYZE ACTIVE SITE CODE ####
##################################

class ACTIVESITE:
    def __init__( self ):
        # get current working directory
        self.working_dir = os.getcwd() + '/'
        
        # make data lists to add over course of program for AS composition
        self.AS_pdb_names = []
        self.AS_lig_res = []
        self.AS_lig_atoms = []
        self.AS_activesite_res = []
        self.AS_activesite_atoms = []
        self.ALA = []
        self.CYS = []
        self.ASP = []
        self.GLU = []
        self.PHE = []
        self.GLY = []
        self.HIS = []
        self.ILE = []
        self.LYS = []
        self.LEU = []
        self.MET = []
        self.ASN = []
        self.PRO = []
        self.GLN = []
        self.ARG = []
        self.SER = []
        self.THR = []
        self.VAL = []
        self.TRP = []
        self.TYR = []
        self.AS_num_activesite_nonpolar_atoms = []
        self.AS_num_activesite_polar_atoms = []
        self.AS_num_activesite_unk_atom_types = []
        self.AS_num_ligand_nonpolar_atoms = []
        self.AS_num_ligand_polar_atoms = []
        self.AS_num_ligand_metal_atoms = []
        self.AS_num_ligand_unk_atom_types = []
        
        # for determining the composition of the total protein
        self.tot_ALA = []
        self.tot_CYS = []
        self.tot_ASP = []
        self.tot_GLU = []
        self.tot_PHE = []
        self.tot_GLY = []
        self.tot_HIS = []
        self.tot_ILE = []
        self.tot_LYS = []
        self.tot_LEU = []
        self.tot_MET = []
        self.tot_ASN = []
        self.tot_PRO = []
        self.tot_GLN = []
        self.tot_ARG = []
        self.tot_SER = []
        self.tot_THR = []
        self.tot_VAL = []
        self.tot_TRP = []
        self.tot_TYR = []
        
        # percentage data holders
        self.percentage_activesite_nonpolar = []
        self.percentage_activesite_polar = []
        self.percentage_ligand_nonpolar = []
        self.percentage_ligand_polar = []
        self.percentage_ligand_metal = []
        self.percentage_activesite_ALA = []
        self.percentage_activesite_CYS = []
        self.percentage_activesite_ASP = []
        self.percentage_activesite_GLU = []
        self.percentage_activesite_PHE = []
        self.percentage_activesite_GLY = []
        self.percentage_activesite_HIS = []
        self.percentage_activesite_ILE = []
        self.percentage_activesite_LYS = []
        self.percentage_activesite_LEU = []
        self.percentage_activesite_MET = []
        self.percentage_activesite_ASN = []
        self.percentage_activesite_PRO = []
        self.percentage_activesite_GLN = []
        self.percentage_activesite_ARG = []
        self.percentage_activesite_SER = []
        self.percentage_activesite_THR = []
        self.percentage_activesite_VAL = []
        self.percentage_activesite_TRP = []
        self.percentage_activesite_TYR = []
        self.percentage_tot_ALA = []
        self.percentage_tot_CYS = []
        self.percentage_tot_ASP = []
        self.percentage_tot_GLU = []
        self.percentage_tot_PHE = []
        self.percentage_tot_GLY = []
        self.percentage_tot_HIS = []
        self.percentage_tot_ILE = []
        self.percentage_tot_LYS = []
        self.percentage_tot_LEU = []
        self.percentage_tot_MET = []
        self.percentage_tot_ASN = []
        self.percentage_tot_PRO = []
        self.percentage_tot_GLN = []
        self.percentage_tot_ARG = []
        self.percentage_tot_SER = []
        self.percentage_tot_THR = []
        self.percentage_tot_VAL = []
        self.percentage_tot_TRP = []
        self.percentage_tot_TYR = []
        
        # make data lists to add over course of program for AA composition per ligand residue in each pdb
        # each pdb name should show up as many times as it has ligand residues that fit the user's criteria
        self.AS_pdb_names_per_lig = []
        self.AS_lig_uniq_res_names_per_lig = []
        self.AS_lig_res_names_per_lig = []
        self.AS_num_lig_atoms_per_lig = []
        self.AS_activesite_res_per_lig = []
        self.AS_activesite_atoms_per_lig = []
        self.ALA_per_lig = []
        self.CYS_per_lig = []
        self.ASP_per_lig = []
        self.GLU_per_lig = []
        self.PHE_per_lig = []
        self.GLY_per_lig = []
        self.HIS_per_lig = []
        self.ILE_per_lig = []
        self.LYS_per_lig = []
        self.LEU_per_lig = []
        self.MET_per_lig = []
        self.ASN_per_lig = []
        self.PRO_per_lig = []
        self.GLN_per_lig = []
        self.ARG_per_lig = []
        self.SER_per_lig = []
        self.THR_per_lig = []
        self.VAL_per_lig = []
        self.TRP_per_lig = []
        self.TYR_per_lig = []
        self.AS_lig_num_ligand_nonpolar_atoms = []
        self.AS_lig_num_ligand_polar_atoms = []
        self.AS_lig_num_ligand_metal_atoms = []
        self.AS_lig_num_ligand_unk_atoms = []
        self.AS_lig_num_activesite_nonpolar_atoms = []
        self.AS_lig_num_activesite_polar_atoms = []
        self.AS_lig_num_activesite_unk_atoms = []

        # percentage data holders
        self.percentage_activesite_per_lig_ALA = []
        self.percentage_activesite_per_lig_CYS = []
        self.percentage_activesite_per_lig_ASP = []
        self.percentage_activesite_per_lig_GLU = []
        self.percentage_activesite_per_lig_PHE = []
        self.percentage_activesite_per_lig_GLY = []
        self.percentage_activesite_per_lig_HIS = []
        self.percentage_activesite_per_lig_ILE = []
        self.percentage_activesite_per_lig_LYS = []
        self.percentage_activesite_per_lig_LEU = []
        self.percentage_activesite_per_lig_MET = []
        self.percentage_activesite_per_lig_ASN = []
        self.percentage_activesite_per_lig_PRO = []
        self.percentage_activesite_per_lig_GLN = []
        self.percentage_activesite_per_lig_ARG = []
        self.percentage_activesite_per_lig_SER = []
        self.percentage_activesite_per_lig_THR = []
        self.percentage_activesite_per_lig_VAL = []
        self.percentage_activesite_per_lig_TRP = []
        self.percentage_activesite_per_lig_TYR = []
        self.percentage_tot_per_lig_ALA = []
        self.percentage_tot_per_lig_CYS = []
        self.percentage_tot_per_lig_ASP = []
        self.percentage_tot_per_lig_GLU = []
        self.percentage_tot_per_lig_PHE = []
        self.percentage_tot_per_lig_GLY = []
        self.percentage_tot_per_lig_HIS = []
        self.percentage_tot_per_lig_ILE = []
        self.percentage_tot_per_lig_LYS = []
        self.percentage_tot_per_lig_LEU = []
        self.percentage_tot_per_lig_MET = []
        self.percentage_tot_per_lig_ASN = []
        self.percentage_tot_per_lig_PRO = []
        self.percentage_tot_per_lig_GLN = []
        self.percentage_tot_per_lig_ARG = []
        self.percentage_tot_per_lig_SER = []
        self.percentage_tot_per_lig_THR = []
        self.percentage_tot_per_lig_VAL = []
        self.percentage_tot_per_lig_TRP = []
        self.percentage_tot_per_lig_TYR = []
        
        
        
    def read_pro_lig_pickles( self, pro_pickle, lig_pickle ):
        # get the four letter code from the passed protein pickle name
        self.name = pro_pickle.split( '/' )[-1][:4].lower()
        
        # collect the corresponding protein and ligand dictionary data
        self.protein = pickle.load( open( pro_pickle, "rb" ) )
        self.ligand = pickle.load( open( lig_pickle, "rb" ) )
        
        # collect number of ligand residues
        self.num_lig_residues = len( self.ligand.keys() )
        
        # collect number of ligand atoms and their atom types
        self.num_lig_atoms = 0
        self.num_lig_nonpolar_atoms = 0
        self.num_lig_polar_atoms = 0
        self.num_lig_metal_atoms = 0
        self.num_lig_unk_atoms = 0
        
        for lig_res in self.ligand.keys():
            self.num_lig_atoms += len( self.ligand[ lig_res ] )
            
            # get the atom lines associated with this residue
            lig_lines = self.ligand[ lig_res ]
            
            # check each element of each atom in the ligand residue
            for line in lig_lines:
                element = line.element
                
                # polar
                if element in polar_atoms:
                    self.num_lig_polar_atoms += 1
                # nonpolar
                elif element in nonpolar_atoms:
                    self.num_lig_nonpolar_atoms += 1
                # metal
                elif element in metal_atoms:
                    self.num_lig_metal_atoms += 1
                # unknown
                else:
                    self.num_lig_unk_atoms += 1
                    
        # collect number of polar, nonpolar, metal, and unk atoms per ligand residue
        self.num_lig_polar_atoms_per_lig = {}
        self.num_lig_nonpolar_atoms_per_lig = {}
        self.num_lig_metal_atoms_per_lig = {}
        self.num_lig_unk_atoms_per_lig = {}
        
        for lig_res in self.ligand.keys():
            # instantiate empty counters
            num_polar = 0
            num_nonpolar = 0
            num_metal = 0
            num_unk = 0
            
            # get the atom lines associated with this residue
            lig_lines = self.ligand[ lig_res ]
            
            # check each element of each atom in the ligand residue
            for line in lig_lines:
                element = line.element
                
                # polar
                if element in polar_atoms:
                    num_polar += 1
                # nonpolar
                elif element in nonpolar_atoms:
                    num_nonpolar += 1
                # metal
                elif element in metal_atoms:
                    num_metal += 1
                # unknown
                else:
                    num_unk += 1
                    
            # add the data to the dictionary for each ligand
            self.num_lig_polar_atoms_per_lig[ lig_res ] = num_polar
            self.num_lig_nonpolar_atoms_per_lig[ lig_res ] = num_nonpolar
            self.num_lig_metal_atoms_per_lig[ lig_res ] = num_metal
            self.num_lig_unk_atoms_per_lig[ lig_res ] = num_unk
            
        return True



    def get_activesite( self, cutoff ):
        # overall activesite dictionary
        # key: unique protein name (resname_reschain_resnum), value: list of ATOM lines per residue
        self.activesite_dict = {}
        
        # unique ligand name (resname_reschain_resnum)
        # value: list of 3-letter amino acid names
        self.activesite_lig_pro_res_name_dict = {}
        
        # unique ligand name (resname_reschain_resnum)
        # value: list of atom_line for each AA ( to get atom count later )
        self.activesite_lig_pro_atoms_dict = {}
        
        # activesite data holders
        self.num_activesite_res = 0
        self.activesite_residues = []
        self.num_activesite_atoms = 0
        self.num_activesite_nonpolar_atoms = 0
        self.num_activesite_polar_atoms = 0
        self.num_activesite_unk_atom_types = 0
        self.activesite_num_nonpolar_atoms = {}
        self.activesite_num_polar_atoms = {}
        self.activesite_num_unk_atoms = {}
        
        for uniq_lig_name in self.ligand.keys():
            # list to store the 3 letter names of all of the protein residues within the cutoff distance of each ligand residue (by unique name)
            AS_names_in_activesite = []
            AS_atoms_in_activesite = []
            
            for hetatm_line in self.ligand[ uniq_lig_name ]:
                # extract coordinates
                x_lig = hetatm_line.x_coord
                y_lig = hetatm_line.y_coord
                z_lig = hetatm_line.z_coord
                lig_xyz = [ x_lig, y_lig, z_lig ]
                
                # for each unique protein residue
                for uniq_pro_name in self.protein.keys():
                    # for each atom in the residue
                    for atom_line in self.protein[ uniq_pro_name ]:
                        # extract coordinates
                        x_pro = atom_line.x_coord
                        y_pro = atom_line.y_coord
                        z_pro = atom_line.z_coord
                        pro_xyz = [ x_pro, y_pro, z_pro ]
                        
                        # check the distance
                        if calc_distance( lig_xyz, pro_xyz ) <= cutoff:
                            # append the line if the unique protein residue has already been counted
                            if uniq_pro_name in self.activesite_residues:
                                if atom_line not in self.activesite_dict[ uniq_pro_name ]:
                                    self.activesite_dict[ uniq_pro_name ].append( atom_line )
                            
                            # store all of the unique names of the protein residues within the activesite
                            # also, store all of the atom_lines for each unique protein residue in the activesite
                            else:
                                self.activesite_residues.append( uniq_pro_name )
                                self.activesite_dict[ uniq_pro_name ] = []
                                self.activesite_dict[ uniq_pro_name ].append( atom_line )
                                
                                # store the 3 letter name of the amino acid within the activesite
                                three_letter_name = uniq_pro_name[ 0:3 ]
                                AS_names_in_activesite.append( three_letter_name )
                                
                            # store atom_line for each unique amino acid within the activesite to get an atom count later
                            if atom_line not in AS_atoms_in_activesite:
                                AS_atoms_in_activesite.append( atom_line )
                        
            # store the list of the 3 letter names for the amino acid within the self.activesite_lig_pro_dict according to which ligand it is near
            self.activesite_lig_pro_res_name_dict[ uniq_lig_name ] = AS_names_in_activesite
            self.activesite_lig_pro_atoms_dict[ uniq_lig_name ] = AS_atoms_in_activesite
            
        # get number of activesite residues
        self.num_activesite_res = len( self.activesite_dict.keys() )
        
        # count the number of activesite atoms
        for uniq_lig_name in self.activesite_lig_pro_atoms_dict.keys():
            self.num_activesite_atoms += len( self.activesite_lig_pro_atoms_dict[ uniq_lig_name ] )
            
            # prepare the counter for nonpolar, polar, and unknown atom types for each unique activesite residue
            self.activesite_num_nonpolar_atoms[ uniq_lig_name ] = 0
            self.activesite_num_polar_atoms[ uniq_lig_name ] = 0
            self.activesite_num_unk_atoms[ uniq_lig_name ] = 0
            
            # count the number of nonpolar and polar atoms
            for pdb_line in self.activesite_lig_pro_atoms_dict[ uniq_lig_name ]:
                # if element is nonpolar
                if pdb_line.element in nonpolar_atoms:
                    # total nonpolar activesite atoms
                    self.num_activesite_nonpolar_atoms += 1
                    
                    # number of nonpolar atoms for this activesite residue
                    self.activesite_num_nonpolar_atoms[ uniq_lig_name ] += 1
                    
                # if element is polar
                elif pdb_line.element in polar_atoms:
                    # total polar activesite atoms
                    self.num_activesite_polar_atoms += 1
                    
                    # number of polar atoms for this activesite residue
                    self.activesite_num_polar_atoms[ uniq_lig_name ] += 1
                    
                # else I don't know what this is
                else:
                    print "      * I didn't know what type of atom", "'%s'" %pdb_line.element(), "is. Please add it to the list"
                    # total unkown activesite atoms
                    self.num_activesite_unk_atom_types += 1
                    
                    # number of unknown atoms for this activesite residue
                    self.activesite_num_unk_atoms[ uniq_lig_name ] += 1
            
        # return information
        # if this ligand has no activesite for some reason
        if len( self.activesite_dict.keys() ) == 0:
            print "## Skipping", self.name, "because it had no activesite residues within", cutoff, "Angstroms of the ligand residue(s) ##"
            return False
        
        # otherwise print relevant information
        else:
            print "  ", self.name, "has", self.num_activesite_res, "activesite residues",
            print "and", self.num_activesite_atoms, "non-hydrogen activesite atoms"
            return True



    def get_activesite_AA_composition( self ):
        # for the activesite composition
        AS_ALA = 0
        AS_CYS = 0
        AS_ASP = 0
        AS_GLU = 0
        AS_PHE = 0
        AS_GLY = 0
        AS_HIS = 0
        AS_ILE = 0
        AS_LYS = 0
        AS_LEU = 0
        AS_MET = 0
        AS_ASN = 0
        AS_PRO = 0
        AS_GLN = 0
        AS_ARG = 0
        AS_SER = 0
        AS_THR = 0
        AS_VAL = 0
        AS_TRP = 0
        AS_TYR = 0

        # for the total protein's composition
        tot_ALA = 0
        tot_CYS = 0
        tot_ASP = 0
        tot_GLU = 0
        tot_PHE = 0
        tot_GLY = 0
        tot_HIS = 0
        tot_ILE = 0
        tot_LYS = 0
        tot_LEU = 0
        tot_MET = 0
        tot_ASN = 0
        tot_PRO = 0
        tot_GLN = 0
        tot_ARG = 0
        tot_SER = 0
        tot_THR = 0
        tot_VAL = 0
        tot_TRP = 0
        tot_TYR = 0

        # loop over each protein residue in the active site and get the first three characters and up the appropriate count
        for pro_res in self.activesite_residues:
            res_name = pro_res[ 0:3 ]
            if res_name == "ALA":
                AS_ALA += 1
            if res_name == "CYS":
                AS_CYS += 1
            if res_name == "ASP":
                AS_ASP += 1
            if res_name == "GLU":
                AS_GLU += 1
            if res_name == "PHE":
                AS_PHE += 1
            if res_name == "GLY":
                AS_GLY += 1
            if res_name == "HIS":
                AS_HIS += 1
            if res_name == "ILE":
                AS_ILE += 1
            if res_name == "LYS":
                AS_LYS += 1
            if res_name == "LEU":
                AS_LEU += 1
            if res_name == "MET":
                AS_MET += 1
            if res_name == "ASN":
                AS_ASN += 1
            if res_name == "PRO":
                AS_PRO += 1
            if res_name == "GLN":
                AS_GLN += 1
            if res_name == "ARG":
                AS_ARG += 1
            if res_name == "SER":
                AS_SER += 1
            if res_name == "THR":
                AS_THR += 1
            if res_name == "VAL":
                AS_VAL += 1
            if res_name == "TRP":
                AS_TRP += 1
            if res_name == "TYR":
                AS_TYR += 1
                
        # loop over each protein residue in the protein and get the first three characters and up the appropriate count
        for pro_res in self.protein.keys():
            res_name = pro_res[ 0:3 ]
            if res_name == "ALA":
                tot_ALA += 1
            if res_name == "CYS":
                tot_CYS += 1
            if res_name == "ASP":
                tot_ASP += 1
            if res_name == "GLU":
                tot_GLU += 1
            if res_name == "PHE":
                tot_PHE += 1
            if res_name == "GLY":
                tot_GLY += 1
            if res_name == "HIS":
                tot_HIS += 1
            if res_name == "ILE":
                tot_ILE += 1
            if res_name == "LYS":
                tot_LYS += 1
            if res_name == "LEU":
                tot_LEU += 1
            if res_name == "MET":
                tot_MET += 1
            if res_name == "ASN":
                tot_ASN += 1
            if res_name == "PRO":
                tot_PRO += 1
            if res_name == "GLN":
                tot_GLN += 1
            if res_name == "ARG":
                tot_ARG += 1
            if res_name == "SER":
                tot_SER += 1
            if res_name == "THR":
                tot_THR += 1
            if res_name == "VAL":
                tot_VAL += 1
            if res_name == "TRP":
                tot_TRP += 1
            if res_name == "TYR":
                tot_TYR += 1
                
        # collect the percentage of polar and nonpolar atoms in the active site
        percentage_activesite_nonpolar = round( 
            float( self.num_activesite_nonpolar_atoms ) / 
            float( self.num_activesite_atoms )
            , 3 )
        percentage_activesite_polar = round( 
            float( self.num_activesite_polar_atoms ) / 
            float( self.num_activesite_atoms )
            , 3 )

        # collect the percentage of polar, nonpolar, and metal atoms in the ligand
        percentage_ligand_nonpolar = round( 
            float( self.num_lig_nonpolar_atoms ) / 
            float( self.num_lig_atoms )
            , 3 )
        percentage_ligand_polar = round( 
            float( self.num_lig_polar_atoms ) / 
            float( self.num_lig_atoms )
            , 3 )
        percentage_ligand_metal = round( 
            float( self.num_lig_metal_atoms ) / 
            float( self.num_lig_atoms )
            , 3 )
        
        # collect the percentage of each specific amino acid in the activesite compared to the total number of amino acids in the activesite
        self.percentage_activesite_ALA.append( round( float( AS_ALA ) / 
                                                      float( len( self.activesite_residues ) )
                                                      , 3 ) )
        self.percentage_activesite_CYS.append( round( float( AS_CYS ) / 
                                                      float( len( self.activesite_residues ) )
                                                      , 3 ) )
        self.percentage_activesite_ASP.append( round( float( AS_ASP ) / 
                                                      float( len( self.activesite_residues ) )
                                                      , 3 ) )
        self.percentage_activesite_GLU.append( round( float( AS_GLU ) / 
                                                      float( len( self.activesite_residues ) )
                                                      , 3 ) )
        self.percentage_activesite_PHE.append( round( float( AS_PHE ) / 
                                                      float( len( self.activesite_residues ) )
                                                      , 3 ) )
        self.percentage_activesite_GLY.append( round( float( AS_GLY ) / 
                                                      float( len( self.activesite_residues ) )
                                                      , 3 ) )
        self.percentage_activesite_HIS.append( round( float( AS_HIS ) / 
                                                      float( len( self.activesite_residues ) )
                                                      , 3 ) )
        self.percentage_activesite_ILE.append( round( float( AS_ILE ) / 
                                                      float( len( self.activesite_residues ) )
                                                      , 3 ) )
        self.percentage_activesite_LYS.append( round( float( AS_LYS ) / 
                                                      float( len( self.activesite_residues ) )
                                                      , 3 ) )
        self.percentage_activesite_LEU.append( round( float( AS_LEU ) / 
                                                      float( len( self.activesite_residues ) )
                                                      , 3 ) )
        self.percentage_activesite_MET.append( round( float( AS_MET ) / 
                                                      float( len( self.activesite_residues ) )
                                                      , 3 ) )
        self.percentage_activesite_ASN.append( round( float( AS_ASN ) / 
                                                      float( len( self.activesite_residues ) )
                                                      , 3 ) )
        self.percentage_activesite_PRO.append( round( float( AS_PRO ) / 
                                                      float( len( self.activesite_residues ) )
                                                      , 3 ) )
        self.percentage_activesite_GLN.append( round( float( AS_GLN ) / 
                                                      float( len( self.activesite_residues ) )
                                                      , 3 ) )
        self.percentage_activesite_ARG.append( round( float( AS_ARG ) / 
                                                      float( len( self.activesite_residues ) )
                                                      , 3 ) )
        self.percentage_activesite_SER.append( round( float( AS_SER ) / 
                                                      float( len( self.activesite_residues ) )
                                                      , 3 ) )
        self.percentage_activesite_THR.append( round( float( AS_THR ) / 
                                                      float( len( self.activesite_residues ) )
                                                      , 3 ) )
        self.percentage_activesite_VAL.append( round( float( AS_VAL ) / 
                                                      float( len( self.activesite_residues ) )
                                                      , 3 ) )
        self.percentage_activesite_TRP.append( round( float( AS_TRP ) / 
                                                      float( len( self.activesite_residues ) )
                                                      , 3 ) )
        self.percentage_activesite_TYR.append( round( float( AS_TYR ) / 
                                                      float( len( self.activesite_residues ) )
                                                      , 3 ) )

        # collect the percentage of each residue in the protein
        self.percentage_tot_ALA.append( round( float( tot_ALA ) / 
                                               float( len( self.protein.keys() ) )
                                               , 3 ) )
        self.percentage_tot_CYS.append( round( float( tot_CYS ) / 
                                               float( len( self.protein.keys() ) )
                                               , 3 ) )
        self.percentage_tot_ASP.append( round( float( tot_ASP ) / 
                                               float( len( self.protein.keys() ) )
                                               , 3 ) )
        self.percentage_tot_GLU.append( round( float( tot_GLU ) / 
                                               float( len( self.protein.keys() ) )
                                               , 3 ) )
        self.percentage_tot_PHE.append( round( float( tot_PHE ) / 
                                               float( len( self.protein.keys() ) )
                                               , 3 ) )
        self.percentage_tot_GLY.append( round( float( tot_GLY ) / 
                                               float( len( self.protein.keys() ) )
                                               , 3 ) )
        self.percentage_tot_HIS.append( round( float( tot_HIS ) / 
                                               float( len( self.protein.keys() ) )
                                               , 3 ) )
        self.percentage_tot_ILE.append( round( float( tot_ILE ) / 
                                               float( len( self.protein.keys() ) )
                                               , 3 ) )
        self.percentage_tot_LYS.append( round( float( tot_LYS ) / 
                                               float( len( self.protein.keys() ) )
                                               , 3 ) )
        self.percentage_tot_LEU.append( round( float( tot_LEU ) / 
                                               float( len( self.protein.keys() ) )
                                               , 3 ) )
        self.percentage_tot_MET.append( round( float( tot_MET ) / 
                                               float( len( self.protein.keys() ) )
                                               , 3 ) )
        self.percentage_tot_ASN.append( round( float( tot_ASN ) / 
                                               float( len( self.protein.keys() ) )
                                               , 3 ) )
        self.percentage_tot_PRO.append( round( float( tot_PRO ) / 
                                               float( len( self.protein.keys() ) )
                                               , 3 ) )
        self.percentage_tot_GLN.append( round( float( tot_GLN ) / 
                                               float( len( self.protein.keys() ) )
                                               , 3 ) )
        self.percentage_tot_ARG.append( round( float( tot_ARG ) / 
                                               float( len( self.protein.keys() ) )
                                               , 3 ) )
        self.percentage_tot_SER.append( round( float( tot_SER ) / 
                                               float( len( self.protein.keys() ) )
                                               , 3 ) )
        self.percentage_tot_THR.append( round( float( tot_THR ) / 
                                               float( len( self.protein.keys() ) )
                                               , 3 ) )
        self.percentage_tot_VAL.append( round( float( tot_VAL ) / 
                                               float( len( self.protein.keys() ) )
                                               , 3 ) )
        self.percentage_tot_TRP.append( round( float( tot_TRP ) / 
                                               float( len( self.protein.keys() ) )
                                               , 3 ) )
        self.percentage_tot_TYR.append( round( float( tot_TYR ) / 
                                               float( len( self.protein.keys() ) )
                                               , 3 ) )

        # append all of the final data to the self.lists
        # because if the analysis got this far, that means there actually is data to collect
        self.AS_pdb_names.append( self.name )
        self.AS_lig_res.append( self.num_lig_residues )
        self.AS_lig_atoms.append( self.num_lig_atoms )
        self.AS_activesite_res.append( self.num_activesite_res )
        self.AS_activesite_atoms.append( self.num_activesite_atoms )
        self.ALA.append( AS_ALA )
        self.CYS.append( AS_CYS )
        self.ASP.append( AS_ASP )
        self.GLU.append( AS_GLU )
        self.PHE.append( AS_PHE )
        self.GLY.append( AS_GLY )
        self.HIS.append( AS_HIS )
        self.ILE.append( AS_ILE )
        self.LYS.append( AS_LYS )
        self.LEU.append( AS_LEU )
        self.MET.append( AS_MET )
        self.ASN.append( AS_ASN )
        self.PRO.append( AS_PRO )
        self.GLN.append( AS_GLN )
        self.ARG.append( AS_ARG )
        self.SER.append( AS_SER )
        self.THR.append( AS_THR )
        self.VAL.append( AS_VAL )
        self.TRP.append( AS_TRP )
        self.TYR.append( AS_TYR )
        self.tot_ALA.append( tot_ALA )
        self.tot_CYS.append( tot_CYS )
        self.tot_ASP.append( tot_ASP )
        self.tot_GLU.append( tot_GLU )
        self.tot_PHE.append( tot_PHE )
        self.tot_GLY.append( tot_GLY )
        self.tot_HIS.append( tot_HIS )
        self.tot_ILE.append( tot_ILE )
        self.tot_LYS.append( tot_LYS )
        self.tot_LEU.append( tot_LEU )
        self.tot_MET.append( tot_MET )
        self.tot_ASN.append( tot_ASN )
        self.tot_PRO.append( tot_PRO )
        self.tot_GLN.append( tot_GLN )
        self.tot_ARG.append( tot_ARG )
        self.tot_SER.append( tot_SER )
        self.tot_THR.append( tot_THR )
        self.tot_VAL.append( tot_VAL )
        self.tot_TRP.append( tot_TRP )
        self.tot_TYR.append( tot_TYR )
        self.AS_num_activesite_nonpolar_atoms.append( self.num_activesite_nonpolar_atoms )
        self.AS_num_activesite_polar_atoms.append( self.num_activesite_polar_atoms )
        self.AS_num_activesite_unk_atom_types.append( self.num_activesite_unk_atom_types )
        self.AS_num_ligand_nonpolar_atoms.append( self.num_lig_nonpolar_atoms )
        self.AS_num_ligand_polar_atoms.append( self.num_lig_polar_atoms )
        self.AS_num_ligand_metal_atoms.append( self.num_lig_metal_atoms )
        self.AS_num_ligand_unk_atom_types.append( self.num_lig_unk_atoms )
        
        return True



    def get_activesite_AA_composition_per_lig_res( self ):
        # goes through each unique ligand residue and counts the number of each amino acid within the cutoff distance around it
        for uniq_lig_name in self.activesite_lig_pro_res_name_dict.keys():
            # append the pdb names to the data list
            self.AS_pdb_names_per_lig.append( self.name )
            
            # append information about each ligand residue
            self.AS_lig_res_names_per_lig.append( uniq_lig_name.split( '_' )[0] )
            self.AS_lig_uniq_res_names_per_lig.append( uniq_lig_name )
            self.AS_num_lig_atoms_per_lig.append( len( self.ligand[ uniq_lig_name ] ) )
            
            # count and append the number of nonpolar and polar ligand atoms
            # this information was collected when reading in the pickle files
            num_nonpolar_lig_atoms = self.num_lig_nonpolar_atoms_per_lig[ uniq_lig_name ]
            num_polar_lig_atoms = self.num_lig_polar_atoms_per_lig[ uniq_lig_name ]
            num_metal_lig_atoms = self.num_lig_metal_atoms_per_lig[ uniq_lig_name ]
            num_unk_lig_atoms = self.num_lig_unk_atoms_per_lig[ uniq_lig_name ]
                    
            self.AS_lig_num_ligand_nonpolar_atoms.append( num_nonpolar_lig_atoms )
            self.AS_lig_num_ligand_polar_atoms.append( num_polar_lig_atoms )
            self.AS_lig_num_ligand_metal_atoms.append( num_metal_lig_atoms )
            self.AS_lig_num_ligand_unk_atoms.append( num_unk_lig_atoms )
            
            # append information about all the activesite residues
            self.AS_activesite_res_per_lig.append( len( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ] ) )
            self.AS_activesite_atoms_per_lig.append( len( self.activesite_lig_pro_atoms_dict[ uniq_lig_name ] ) )
            
            # count and append the number of nonpolar and polar activesite atoms
            num_nonpolar_activesite_atoms = self.activesite_num_nonpolar_atoms[ uniq_lig_name ]
            num_polar_activesite_atoms = self.activesite_num_polar_atoms[ uniq_lig_name ]
            num_unk_activesite_atoms = self.activesite_num_unk_atoms[ uniq_lig_name ]
            
            self.AS_lig_num_activesite_nonpolar_atoms.append( num_nonpolar_activesite_atoms )
            self.AS_lig_num_activesite_polar_atoms.append( num_polar_activesite_atoms )
            self.AS_lig_num_activesite_unk_atoms.append( num_unk_activesite_atoms )
            
            # count the number of amino acid residues around the ligand and append to data lists
            self.ALA_per_lig.append( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "ALA" ) )
            self.CYS_per_lig.append( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "CYS" ) )
            self.ASP_per_lig.append( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "ASP" ) )
            self.GLU_per_lig.append( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "GLU" ) )
            self.PHE_per_lig.append( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "PHE" ) )
            self.GLY_per_lig.append( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "GLY" ) )
            self.HIS_per_lig.append( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "HIS" ) )
            self.ILE_per_lig.append( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "ILE" ) )
            self.LYS_per_lig.append( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "LYS" ) )
            self.LEU_per_lig.append( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "LEU" ) )
            self.MET_per_lig.append( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "MET" ) )
            self.ASN_per_lig.append( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "ASN" ) )
            self.PRO_per_lig.append( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "PRO" ) )
            self.GLN_per_lig.append( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "GLN" ) )
            self.ARG_per_lig.append( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "ARG" ) )
            self.SER_per_lig.append( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "SER" ) )
            self.THR_per_lig.append( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "THR" ) )
            self.VAL_per_lig.append( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "VAL" ) )
            self.TRP_per_lig.append( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "TRP" ) )
            self.TYR_per_lig.append( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "TYR" ) )
            
            # collect percentage data of number of a specific amino acid in the ligand's activesite versus the number of total amino acids in that activesite
            self.percentage_activesite_per_lig_ALA.append( 
                round( 
                    float( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "ALA" ) ) / 
                    float( len( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ] ) )
                    , 3 ) )
            self.percentage_activesite_per_lig_CYS.append( 
                round( 
                    float( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "CYS" ) ) / 
                    float( len( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ] ) )
                    , 3 ) )
            self.percentage_activesite_per_lig_ASP.append( 
                round( 
                    float( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "ASP" ) ) / 
                    float( len( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ] ) )
                    , 3 ) )
            self.percentage_activesite_per_lig_GLU.append( 
                round( 
                    float( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "GLU" ) ) / 
                    float( len( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ] ) )
                    , 3 ) )
            self.percentage_activesite_per_lig_PHE.append( 
                round( 
                    float( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "PHE" ) ) / 
                    float( len( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ] ) )
                    , 3 ) )
            self.percentage_activesite_per_lig_GLY.append( 
                round( 
                    float( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "GLY" ) ) / 
                    float( len( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ] ) )
                    , 3 ) )
            self.percentage_activesite_per_lig_HIS.append( 
                round( 
                    float( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "HIS" ) ) / 
                    float( len( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ] ) )
                    , 3 ) )
            self.percentage_activesite_per_lig_ILE.append( 
                round( 
                    float( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "ILE" ) ) / 
                    float( len( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ] ) )
                    , 3 ) )
            self.percentage_activesite_per_lig_LYS.append( 
                round( 
                    float( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "LYS" ) ) / 
                    float( len( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ] ) )
                    , 3 ) )
            self.percentage_activesite_per_lig_LEU.append( 
                round( 
                    float( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "LEU" ) ) / 
                    float( len( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ] ) )
                    , 3 ) )
            self.percentage_activesite_per_lig_MET.append( 
                round( 
                    float( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "MET" ) ) / 
                    float( len( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ] ) )
                    , 3 ) )
            self.percentage_activesite_per_lig_ASN.append( 
                round( 
                    float( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "ASN" ) ) / 
                    float( len( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ] ) )
                    , 3 ) )
            self.percentage_activesite_per_lig_PRO.append( 
                round( 
                    float( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "PRO" ) ) / 
                    float( len( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ] ) )
                    , 3 ) )
            self.percentage_activesite_per_lig_GLN.append( 
                round( 
                    float( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "GLN" ) ) / 
                    float( len( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ] ) )
                    , 3 ) )
            self.percentage_activesite_per_lig_ARG.append( 
                round( 
                    float( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "ARG" ) ) / 
                    float( len( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ] ) )
                    , 3 ) )
            self.percentage_activesite_per_lig_SER.append( 
                round( 
                    float( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "SER" ) ) / 
                    float( len( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ] ) )
                    , 3 ) )
            self.percentage_activesite_per_lig_THR.append( 
                round( 
                    float( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "THR" ) ) / 
                    float( len( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ] ) )
                    , 3 ) )
            self.percentage_activesite_per_lig_VAL.append( 
                round( 
                    float( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "VAL" ) ) / 
                    float( len( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ] ) )
                    , 3 ) )
            self.percentage_activesite_per_lig_TRP.append( 
                round( 
                    float( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "TRP" ) ) / 
                    float( len( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ] ) )
                    , 3 ) )
            self.percentage_activesite_per_lig_TYR.append( 
                round( 
                    float( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ].count( "TYR" ) ) / 
                    float( len( self.activesite_lig_pro_res_name_dict[ uniq_lig_name ] ) )
                    , 3 ) )
            
            # recollect the percentage of each residue in the protein
            self.percentage_tot_per_lig_ALA.append( self.percentage_tot_ALA[-1] )
            self.percentage_tot_per_lig_CYS.append( self.percentage_tot_CYS[-1] )
            self.percentage_tot_per_lig_ASP.append( self.percentage_tot_ASP[-1] )
            self.percentage_tot_per_lig_GLU.append( self.percentage_tot_GLU[-1] )
            self.percentage_tot_per_lig_PHE.append( self.percentage_tot_PHE[-1] )
            self.percentage_tot_per_lig_GLY.append( self.percentage_tot_GLY[-1] )
            self.percentage_tot_per_lig_HIS.append( self.percentage_tot_HIS[-1] )
            self.percentage_tot_per_lig_ILE.append( self.percentage_tot_ILE[-1] )
            self.percentage_tot_per_lig_LYS.append( self.percentage_tot_LYS[-1] )
            self.percentage_tot_per_lig_LEU.append( self.percentage_tot_LEU[-1] )
            self.percentage_tot_per_lig_MET.append( self.percentage_tot_MET[-1] )
            self.percentage_tot_per_lig_ASN.append( self.percentage_tot_ASN[-1] )
            self.percentage_tot_per_lig_PRO.append( self.percentage_tot_PRO[-1] )
            self.percentage_tot_per_lig_GLN.append( self.percentage_tot_GLN[-1] )
            self.percentage_tot_per_lig_ARG.append( self.percentage_tot_ARG[-1] )
            self.percentage_tot_per_lig_SER.append( self.percentage_tot_SER[-1] )
            self.percentage_tot_per_lig_THR.append( self.percentage_tot_THR[-1] )
            self.percentage_tot_per_lig_VAL.append( self.percentage_tot_VAL[-1] )
            self.percentage_tot_per_lig_TRP.append( self.percentage_tot_TRP[-1] )
            self.percentage_tot_per_lig_TYR.append( self.percentage_tot_TYR[-1] )
                                    
        return True
