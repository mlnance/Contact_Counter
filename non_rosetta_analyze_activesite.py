#!/usr/bin/python
__author__ = "morganlnance"


#####
# path to pymol executable
# must change for each machine!
#pymol_exe_dir = "/Applications/MacPyMOL.app/Contents/MacOS/MacPyMOL"
#####


'''
TODO: finish fixing both of the per_lig functions. Should just append data to the lists directly after counting the stuff because otherwise it's just too complicated to keep track of...maybe should fix the other two functions too to do this. Ie. get rid of collect_activesite_data
'''


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
        
        # make data lists to add over course of program for activesite (AS) composition
        self.AS_pdb_names_list = [] 
        self.AS_lig_res_list = []
        self.AS_lig_atoms_list = []
        self.AS_activesite_res_list = []
        self.AS_activesite_atoms_list = []
        self.AS_num_activesite_nonpolar_atoms_types_list = []
        self.AS_num_activesite_polar_atom_types_list = []
        self.AS_num_activesite_unk_atom_types_list = []
        self.AS_num_ligand_nonpolar_atom_types_list = []
        self.AS_num_ligand_polar_atom_types_list = []
        self.AS_num_ligand_metal_atom_types_list = []
        self.AS_num_ligand_unk_atom_types_list = []
        
        # for the specific AA composition in the activesite
        self.AS_ALA_list = []
        self.AS_CYS_list = []
        self.AS_ASP_list = []
        self.AS_GLU_list = []
        self.AS_PHE_list = []
        self.AS_GLY_list = []
        self.AS_HIS_list = []
        self.AS_ILE_list = []
        self.AS_LYS_list = []
        self.AS_LEU_list = []
        self.AS_MET_list = []
        self.AS_ASN_list = []
        self.AS_PRO_list = []
        self.AS_GLN_list = []
        self.AS_ARG_list = []
        self.AS_SER_list = []
        self.AS_THR_list = []
        self.AS_VAL_list = []
        self.AS_TRP_list = []
        self.AS_TYR_list = []
        
        # for determining the composition of the total protein
        self.tot_ALA_list = []
        self.tot_CYS_list = []
        self.tot_ASP_list = []
        self.tot_GLU_list = []
        self.tot_PHE_list = []
        self.tot_GLY_list = []
        self.tot_HIS_list = []
        self.tot_ILE_list = []
        self.tot_LYS_list = []
        self.tot_LEU_list = []
        self.tot_MET_list = []
        self.tot_ASN_list = []
        self.tot_PRO_list = []
        self.tot_GLN_list = []
        self.tot_ARG_list = []
        self.tot_SER_list = []
        self.tot_THR_list = []
        self.tot_VAL_list = []
        self.tot_TRP_list = []
        self.tot_TYR_list = []
        
        # percentage data holders
        self.percentage_activesite_nonpolar_list = []
        self.percentage_activesite_polar_list = []
        self.percentage_ligand_nonpolar_list = []
        self.percentage_ligand_polar_list = []
        self.percentage_ligand_metal_list = []
        self.percentage_activesite_ALA_list = []
        self.percentage_activesite_CYS_list = []
        self.percentage_activesite_ASP_list = []
        self.percentage_activesite_GLU_list = []
        self.percentage_activesite_PHE_list = []
        self.percentage_activesite_GLY_list = []
        self.percentage_activesite_HIS_list = []
        self.percentage_activesite_ILE_list = []
        self.percentage_activesite_LYS_list = []
        self.percentage_activesite_LEU_list = []
        self.percentage_activesite_MET_list = []
        self.percentage_activesite_ASN_list = []
        self.percentage_activesite_PRO_list = []
        self.percentage_activesite_GLN_list = []
        self.percentage_activesite_ARG_list = []
        self.percentage_activesite_SER_list = []
        self.percentage_activesite_THR_list = []
        self.percentage_activesite_VAL_list = []
        self.percentage_activesite_TRP_list = []
        self.percentage_activesite_TYR_list = []
        self.percentage_tot_ALA_list = []
        self.percentage_tot_CYS_list = []
        self.percentage_tot_ASP_list = []
        self.percentage_tot_GLU_list = []
        self.percentage_tot_PHE_list = []
        self.percentage_tot_GLY_list = []
        self.percentage_tot_HIS_list = []
        self.percentage_tot_ILE_list = []
        self.percentage_tot_LYS_list = []
        self.percentage_tot_LEU_list = []
        self.percentage_tot_MET_list = []
        self.percentage_tot_ASN_list = []
        self.percentage_tot_PRO_list = []
        self.percentage_tot_GLN_list = []
        self.percentage_tot_ARG_list = []
        self.percentage_tot_SER_list = []
        self.percentage_tot_THR_list = []
        self.percentage_tot_VAL_list = []
        self.percentage_tot_TRP_list = []
        self.percentage_tot_TYR_list = []
        
        
        # make data lists to add over course of program for AA composition per ligand residue in each PDB
        # each PDB name should show up as many times as it has ligand residues that fit the user's criteria
        self.AS_pdb_names_per_lig_list = []
        self.AS_lig_uniq_res_names_per_lig_list = []
        self.AS_lig_res_names_per_lig_list = []
        self.AS_num_lig_atoms_per_lig_list = []
        self.AS_activesite_res_per_lig_list = []
        self.AS_activesite_atoms_per_lig_list = []
        self.AS_lig_num_ligand_nonpolar_atoms_per_lig_list = []
        self.AS_lig_num_ligand_polar_atoms_per_lig_list = []
        self.AS_lig_num_ligand_metal_atoms_per_lig_list = []
        self.AS_lig_num_ligand_unk_atoms_per_lig_list = []
        self.AS_lig_num_activesite_nonpolar_atoms_per_lig_list = []
        self.AS_lig_num_activesite_polar_atoms_per_lig_list = []
        self.AS_lig_num_activesite_unk_atoms_per_lig_list = []
        
        # for the specific AA composition in each activesite
        self.AS_ALA_per_lig_list = []
        self.AS_CYS_per_lig_list = []
        self.AS_ASP_per_lig_list = []
        self.AS_GLU_per_lig_list = []
        self.AS_PHE_per_lig_list = []
        self.AS_GLY_per_lig_list = []
        self.AS_HIS_per_lig_list = []
        self.AS_ILE_per_lig_list = []
        self.AS_LYS_per_lig_list = []
        self.AS_LEU_per_lig_list = []
        self.AS_MET_per_lig_list = []
        self.AS_ASN_per_lig_list = []
        self.AS_PRO_per_lig_list = []
        self.AS_GLN_per_lig_list = []
        self.AS_ARG_per_lig_list = []
        self.AS_SER_per_lig_list = []
        self.AS_THR_per_lig_list = []
        self.AS_VAL_per_lig_list = []
        self.AS_TRP_per_lig_list = []
        self.AS_TYR_per_lig_list = []

        # percentage data holders per ligand
        self.percentage_activesite_per_lig_ALA_list = []
        self.percentage_activesite_per_lig_CYS_list = []
        self.percentage_activesite_per_lig_ASP_list = []
        self.percentage_activesite_per_lig_GLU_list = []
        self.percentage_activesite_per_lig_PHE_list = []
        self.percentage_activesite_per_lig_GLY_list = []
        self.percentage_activesite_per_lig_HIS_list = []
        self.percentage_activesite_per_lig_ILE_list = []
        self.percentage_activesite_per_lig_LYS_list = []
        self.percentage_activesite_per_lig_LEU_list = []
        self.percentage_activesite_per_lig_MET_list = []
        self.percentage_activesite_per_lig_ASN_list = []
        self.percentage_activesite_per_lig_PRO_list = []
        self.percentage_activesite_per_lig_GLN_list = []
        self.percentage_activesite_per_lig_ARG_list = []
        self.percentage_activesite_per_lig_SER_list = []
        self.percentage_activesite_per_lig_THR_list = []
        self.percentage_activesite_per_lig_VAL_list = []
        self.percentage_activesite_per_lig_TRP_list = []
        self.percentage_activesite_per_lig_TYR_list = []
        self.percentage_tot_per_lig_ALA_list = []
        self.percentage_tot_per_lig_CYS_list = []
        self.percentage_tot_per_lig_ASP_list = []
        self.percentage_tot_per_lig_GLU_list = []
        self.percentage_tot_per_lig_PHE_list = []
        self.percentage_tot_per_lig_GLY_list = []
        self.percentage_tot_per_lig_HIS_list = []
        self.percentage_tot_per_lig_ILE_list = []
        self.percentage_tot_per_lig_LYS_list = []
        self.percentage_tot_per_lig_LEU_list = []
        self.percentage_tot_per_lig_MET_list = []
        self.percentage_tot_per_lig_ASN_list = []
        self.percentage_tot_per_lig_PRO_list = []
        self.percentage_tot_per_lig_GLN_list = []
        self.percentage_tot_per_lig_ARG_list = []
        self.percentage_tot_per_lig_SER_list = []
        self.percentage_tot_per_lig_THR_list = []
        self.percentage_tot_per_lig_VAL_list = []
        self.percentage_tot_per_lig_TRP_list = []
        self.percentage_tot_per_lig_TYR_list = []
        
        
        
    def read_pro_lig_pickles( self, pro_pickle, lig_pickle ):
        # get the four letter code from the passed protein pickle name
        self.name = pro_pickle.split( '/' )[-1][:4].lower()
        self.AS_pdb_names_list.append( self.name )
        
        # collect the corresponding protein and ligand dictionary data
        self.protein = pickle.load( open( pro_pickle, "rb" ) )
        self.ligand = pickle.load( open( lig_pickle, "rb" ) )
        
        return True
    
    
    
    def analyze_ligand( self ):
        # instantiate data holders
        self.num_lig_residues = 0
        self.num_lig_atoms = 0
        self.num_lig_nonpolar_atoms = 0
        self.num_lig_polar_atoms = 0
        self.num_lig_metal_atoms = 0
        self.num_lig_unk_atoms = 0
        self.num_lig_atoms_per_lig = {}
        self.num_lig_polar_atoms_per_lig = {}
        self.num_lig_nonpolar_atoms_per_lig = {}
        self.num_lig_metal_atoms_per_lig = {}
        self.num_lig_unk_atoms_per_lig = {}
                
        # collect number of ligands
        self.num_lig_residues = len( self.ligand.keys() )
        
        # collect overall number of ligand atoms and their atom types
        for lig_res in self.ligand.keys():
            # add to total number of ligand atoms
            self.num_lig_atoms += len( self.ligand[ lig_res ] )
            
            # get the atom lines associated with this residue to count atom types
            lig_lines = self.ligand[ lig_res ]
            
            # check each element of each atom in the ligand
            for lig_line in lig_lines:
                # if polar
                if lig_line.element in polar_atoms:
                    self.num_lig_polar_atoms += 1
                # if nonpolar
                elif lig_line.element in nonpolar_atoms:
                    self.num_lig_nonpolar_atoms += 1
                # if metal
                elif lig_lin.element in metal_atoms:
                    self.num_lig_metal_atoms += 1
                # else unknown
                else:
                    self.num_lig_unk_atoms += 1
                    
                    
        # collect number of polar, nonpolar, metal, and unk atoms of each ligand
        for lig_res in self.ligand.keys():
            # instantiate empty counters for this ligand
            num_polar = 0
            num_nonpolar = 0
            num_metal = 0
            num_unk = 0
            
            # get the atom lines associated with this ligand
            lig_lines = self.ligand[ lig_res ]
            
            # number of atoms in total of this ligand
            num_lig_atoms = len( lig_lines )
            
            # check each element of each atom in the ligand
            for lig_line in lig_lines:
                # else polar
                if lig_line.element in polar_atoms:
                    num_polar += 1
                # else nonpolar
                elif lig_line.element in nonpolar_atoms:
                    num_nonpolar += 1
                # else metal
                elif lig_line.element in metal_atoms:
                    num_metal += 1
                # else unknown
                else:
                    num_unk += 1
                    
            # add the data to the dictionary for each ligand
            self.num_lig_atoms_per_lig[ lig_res ] = num_lig_atoms
            self.num_lig_polar_atoms_per_lig[ lig_res ] = num_polar
            self.num_lig_nonpolar_atoms_per_lig[ lig_res ] = num_nonpolar
            self.num_lig_metal_atoms_per_lig[ lig_res ] = num_metal
            self.num_lig_unk_atoms_per_lig[ lig_res ] = num_unk
            
        # collect the percentage of polar, nonpolar, and metal atoms in the ligand
        self.percentage_ligand_nonpolar = round( 
            float( self.num_lig_nonpolar_atoms ) / 
            float( self.num_lig_atoms )
            , 3 )
        self.percentage_ligand_polar = round( 
            float( self.num_lig_polar_atoms ) / 
            float( self.num_lig_atoms )
            , 3 )
        self.percentage_ligand_metal = round( 
            float( self.num_lig_metal_atoms ) / 
            float( self.num_lig_atoms )
            , 3 )
        
        return True
    
    
    
    def analyze_protein( self ):
        # for the size of the protein
        self.size_of_protein = len( self.protein.keys() )
        
        # for the total protein's composition
        self.tot_ALA = 0
        self.tot_CYS = 0
        self.tot_ASP = 0
        self.tot_GLU = 0
        self.tot_PHE = 0
        self.tot_GLY = 0
        self.tot_HIS = 0
        self.tot_ILE = 0
        self.tot_LYS = 0
        self.tot_LEU = 0
        self.tot_MET = 0
        self.tot_ASN = 0
        self.tot_PRO = 0
        self.tot_GLN = 0
        self.tot_ARG = 0
        self.tot_SER = 0
        self.tot_THR = 0
        self.tot_VAL = 0
        self.tot_TRP = 0
        self.tot_TYR = 0
        
        # for each unique residue in the protein
        for pro_res in self.protein.keys():
            # resname_reschain_reseqpos
            # so splitting on '_' will give the residue name in the first element
            res_name = pro_res.split( '_' )[0]
            
            # count each residue
            if res_name == "ALA":
                self.tot_ALA += 1
            if res_name == "CYS":
                self.tot_CYS += 1
            if res_name == "ASP":
                self.tot_ASP += 1
            if res_name == "GLU":
                self.tot_GLU += 1
            if res_name == "PHE":
                self.tot_PHE += 1
            if res_name == "GLY":
                self.tot_GLY += 1
            if res_name == "HIS":
                self.tot_HIS += 1
            if res_name == "ILE":
                self.tot_ILE += 1
            if res_name == "LYS":
                self.tot_LYS += 1
            if res_name == "LEU":
                self.tot_LEU += 1
            if res_name == "MET":
                self.tot_MET += 1
            if res_name == "ASN":
                self.tot_ASN += 1
            if res_name == "PRO":
                self.tot_PRO += 1
            if res_name == "GLN":
                self.tot_GLN += 1
            if res_name == "ARG":
                self.tot_ARG += 1
            if res_name == "SER":
                self.tot_SER += 1
            if res_name == "THR":
                self.tot_THR += 1
            if res_name == "VAL":
                self.tot_VAL += 1
            if res_name == "TRP":
                self.tot_TRP += 1
            if res_name == "TYR":
                self.tot_TYR += 1
                
        # collect the percentage of each residue in the overall protein
        self.percentage_tot_ALA = ( round( float( self.tot_ALA ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
        self.percentage_tot_CYS = ( round( float( self.tot_CYS ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
        self.percentage_tot_ASP = ( round( float( self.tot_ASP ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
        self.percentage_tot_GLU = ( round( float( self.tot_GLU ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
        self.percentage_tot_PHE = ( round( float( self.tot_PHE ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
        self.percentage_tot_GLY = ( round( float( self.tot_GLY ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
        self.percentage_tot_HIS = ( round( float( self.tot_HIS ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
        self.percentage_tot_ILE = ( round( float( self.tot_ILE ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
        self.percentage_tot_LYS = ( round( float( self.tot_LYS ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
        self.percentage_tot_LEU = ( round( float( self.tot_LEU ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
        self.percentage_tot_MET = ( round( float( self.tot_MET ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
        self.percentage_tot_ASN = ( round( float( self.tot_ASN ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
        self.percentage_tot_PRO = ( round( float( self.tot_PRO ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
        self.percentage_tot_GLN = ( round( float( self.tot_GLN ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
        self.percentage_tot_ARG = ( round( float( self.tot_ARG ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
        self.percentage_tot_SER = ( round( float( self.tot_SER ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
        self.percentage_tot_THR = ( round( float( self.tot_THR ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
        self.percentage_tot_VAL = ( round( float( self.tot_VAL ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
        self.percentage_tot_TRP = ( round( float( self.tot_TRP ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
        self.percentage_tot_TYR = ( round( float( self.tot_TYR ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
        
        return True
    
    
    
    def get_activesite( self, cutoff ):
        # key: unique ligand name (resname_reschain_resnum)
        # value: unique protein names or protein atom lines
        self.activesite_residues_per_lig = {}
        self.activesite_atoms_per_lig = {}
        
        # list of unique activesite residues and atoms
        self.uniq_activesite_residues = []
        self.uniq_activesite_atoms = []
        
        # activesite data holders
        self.num_activesite_residues = 0
        self.num_activesite_atoms = 0
        
        
        # for each unique ligand residue
        for uniq_lig_name in self.ligand.keys():
            # add the unique ligand residue name to self.activesite_dict
            self.activesite_residues_per_lig[ uniq_lig_name ] = []
            self.activesite_atoms_per_lig[ uniq_lig_name ] = []
                
            # for each atom in the ligand residue
            for hetatm_line in self.ligand[ uniq_lig_name ]:
                # extract coordinates
                x_lig = hetatm_line.x_coord
                y_lig = hetatm_line.y_coord
                z_lig = hetatm_line.z_coord
                lig_xyz = [ x_lig, y_lig, z_lig ]
                
                # for each unique protein residue
                for uniq_pro_name in self.protein.keys():
                    # for each atom in the protein residue
                    for atom_line in self.protein[ uniq_pro_name ]:
                        # extract coordinates
                        x_pro = atom_line.x_coord
                        y_pro = atom_line.y_coord
                        z_pro = atom_line.z_coord
                        pro_xyz = [ x_pro, y_pro, z_pro ]
                        
                        # check the distance
                        if calc_distance( lig_xyz, pro_xyz ) <= cutoff:
                            # add the unique protein name to the activesite_residues dictionary
                            if uniq_pro_name not in self.activesite_residues_per_lig[ uniq_lig_name ]:
                                self.activesite_residues_per_lig[ uniq_lig_name ].append( uniq_pro_name )
                                
                            # add the unique protein name to the overall activesite residue list
                            if uniq_pro_name not in self.uniq_activesite_residues:
                                self.uniq_activesite_residues.append( uniq_pro_name )
                                
                            # add the protein atom line to the activesite_atoms dictionary
                            if atom_line not in self.activesite_atoms_per_lig[ uniq_lig_name ]:
                                self.activesite_atoms_per_lig[ uniq_lig_name ].append( atom_line )

                            # add the unique protein name to the overall activesite residue list
                            if atom_line not in self.uniq_activesite_atoms:
                                self.uniq_activesite_atoms.append( atom_line )
                        

        # get number of activesite residues
        self.num_activesite_residues = len( self.uniq_activesite_residues )
        
        # get number of activesite atoms
        # this is the count of atoms actually in the activesite
        # not just the atoms of the residues in the activesite
        self.num_activesite_atoms = len( self.uniq_activesite_atoms )
        
        # return information
        # if this ligand has no activesite for some reason
        if self.num_activesite_residues == 0:
            print "## Skipping", self.name, "because it had no activesite residues within", cutoff, "Angstroms of the ligand residue(s) ##"
            
            return False
        
        # otherwise print relevant information
        else:
            print "  ", self.name, "has", self.num_activesite_residues, "activesite residues",
            print "and", self.num_activesite_atoms, "non-hydrogen activesite atoms"
            
            return True



    def get_activesite_atom_composition( self ):
        # instantiate data holders
        self.num_activesite_nonpolar_atom_types = 0
        self.num_activesite_polar_atom_types = 0
        self.num_activesite_unk_atom_types = 0
        self.num_activesite_nonpolar_atoms = {}
        self.num_activesite_polar_atoms = {}
        self.num_activesite_unk_atoms = {}
        
        # collect the types of the unique atoms in the activesite
        for atom_line in self.uniq_activesite_atoms:
            # if nonpolar
            if atom_line.element in nonpolar_atoms:
                self.num_activesite_nonpolar_atom_types += 1
            # if polar
            elif atom_line.element in polar_atoms:
                self.num_activesite_polar_atom_types += 1
            # else unknown
            else:
                print "      * I didn't know what type of atom", "'%s'" %atom_line.element, "is. Please add it to the list"
                self.num_activesite_unk_atom_types += 1
                
        # add the data to the dictionaries given the PDB name
        self.num_activesite_nonpolar_atoms[ self.name ] =self.num_activesite_nonpolar_atom_types
        self.num_activesite_polar_atoms[ self.name ] =self.num_activesite_polar_atom_types
        self.num_activesite_unk_atoms[ self.name ] =self.num_activesite_unk_atom_types

        # collect the percentage of polar and nonpolar atoms in the activesite
        self.percentage_activesite_nonpolar = round( 
            float( self.num_activesite_nonpolar_atom_types ) / 
            float( self.num_activesite_atoms )
            , 3 )
        self.percentage_activesite_polar = round( 
            float( self.num_activesite_polar_atom_types ) / 
            float( self.num_activesite_atoms )
            , 3 )
                
        return True
    
    
    
    def get_activesite_atom_composition_per_lig( self ):
        # for each ligand residue, collect the number of each atom type
        for uniq_lig_res in self.activesite_atoms_per_lig.keys():
            # add some PDB and ligand data to the data lists
            self.AS_pdb_names_per_lig_list.append( self.name )
            self.AS_lig_uniq_res_names_per_lig_list.append( uniq_lig_res )
            self.AS_lig_res_names_per_lig_list.append( uniq_lig_res.split( '_' )[0] )
            self.AS_num_lig_atoms_per_lig_list.append( 
                self.num_lig_atoms_per_lig[ uniq_lig_res ] )
            self.AS_lig_num_ligand_nonpolar_atoms_per_lig_list.append( 
                self.num_lig_nonpolar_atoms_per_lig[ uniq_lig_res ] )
            self.AS_lig_num_ligand_polar_atoms_per_lig_list.append( 
                self.num_lig_polar_atoms_per_lig[ uniq_lig_res ] )
            self.AS_lig_num_ligand_metal_atoms_per_lig_list.append( 
                self.num_lig_metal_atoms_per_lig[ uniq_lig_res ] )
            self.AS_lig_num_ligand_unk_atoms_per_lig_list.append( 
                self.num_lig_unk_atoms_per_lig[ uniq_lig_res ] )

            # instantiate empty counters for each ligand
            num_activesite_nonpolar_atom_types_per_lig = 0
            num_activesite_polar_atom_types_per_lig = 0
            num_activesite_unk_atom_types_per_lig = 0
            
            # check the element of each atom in the activesite around this ligand
            activesite_atom_lines = self.activesite_atoms_per_lig[ uniq_lig_res ]
            num_activesite_atoms_per_lig = len( activesite_atom_lines )
            
            for atom_line in activesite_atom_lines:
                # if nonpolar
                if atom_line.element in nonpolar_atoms:
                    num_activesite_nonpolar_atom_types_per_lig += 1
                # if polar
                elif atom_line.element in polar_atoms:
                    num_activesite_polar_atom_types_per_lig += 1
                # else unknown
                else:
                    print "      * I didn't know what type of atom", "'%s'" %atom_line.element, "is. Please add it to the list"
                    num_activesite_unk_atom_types_per_lig += 1
                    
            # add the data to the dictionaries for each ligand residue
            self.AS_activesite_atoms_per_lig_list.append( num_activesite_atoms_per_lig )
            self.AS_lig_num_activesite_nonpolar_atoms_per_lig_list.append( num_activesite_nonpolar_atom_types_per_lig )
            self.AS_lig_num_activesite_polar_atoms_per_lig_list.append( num_activesite_polar_atom_types_per_lig )
            self.AS_lig_num_activesite_unk_atoms_per_lig_list.append( num_activesite_unk_atom_types_per_lig )
            
        return True
    
    
    
    def get_activesite_AA_composition( self ):
        # for the overall activesite composition
        self.AS_ALA = 0
        self.AS_CYS = 0
        self.AS_ASP = 0
        self.AS_GLU = 0
        self.AS_PHE = 0
        self.AS_GLY = 0
        self.AS_HIS = 0
        self.AS_ILE = 0
        self.AS_LYS = 0
        self.AS_LEU = 0
        self.AS_MET = 0
        self.AS_ASN = 0
        self.AS_PRO = 0
        self.AS_GLN = 0
        self.AS_ARG = 0
        self.AS_SER = 0
        self.AS_THR = 0
        self.AS_VAL = 0
        self.AS_TRP = 0
        self.AS_TYR = 0

        # for each unique residue in the activesite
        for pro_res in self.uniq_activesite_residues:
            # resname_reschain_reseqpos
            # so splitting on '_' will give the residue name in the first element
            res_name = pro_res.split( '_' )[0]
            
            # count the residue appearances in all of the activesites
            if res_name == "ALA":
                self.AS_ALA += 1
            if res_name == "CYS":
                self.AS_CYS += 1
            if res_name == "ASP":
                self.AS_ASP += 1
            if res_name == "GLU":
                self.AS_GLU += 1
            if res_name == "PHE":
                self.AS_PHE += 1
            if res_name == "GLY":
                self.AS_GLY += 1
            if res_name == "HIS":
                self.AS_HIS += 1
            if res_name == "ILE":
                self.AS_ILE += 1
            if res_name == "LYS":
                self.AS_LYS += 1
            if res_name == "LEU":
                self.AS_LEU += 1
            if res_name == "MET":
                self.AS_MET += 1
            if res_name == "ASN":
                self.AS_ASN += 1
            if res_name == "PRO":
                self.AS_PRO += 1
            if res_name == "GLN":
                self.AS_GLN += 1
            if res_name == "ARG":
                self.AS_ARG += 1
            if res_name == "SER":
                self.AS_SER += 1
            if res_name == "THR":
                self.AS_THR += 1
            if res_name == "VAL":
                self.AS_VAL += 1
            if res_name == "TRP":
                self.AS_TRP += 1
            if res_name == "TYR":
                self.AS_TYR += 1
        
        # collect the percentage of each specific amino acid in the activesite compared to the total number of amino acids in the activesite
        self.percentage_activesite_ALA = ( round( float( self.AS_ALA ) / 
                                                  float( self.num_activesite_residues )
                                                  , 3 ) )
        self.percentage_activesite_CYS = ( round( float( self.AS_CYS ) / 
                                                  float( self.num_activesite_residues )
                                                  , 3 ) )
        self.percentage_activesite_ASP = ( round( float( self.AS_ASP ) / 
                                                  float( self.num_activesite_residues )
                                                  , 3 ) )
        self.percentage_activesite_GLU = ( round( float( self.AS_GLU ) / 
                                                  float( self.num_activesite_residues )
                                                  , 3 ) )
        self.percentage_activesite_PHE = ( round( float( self.AS_PHE ) / 
                                                  float( self.num_activesite_residues )
                                                  , 3 ) )
        self.percentage_activesite_GLY = ( round( float( self.AS_GLY ) / 
                                                  float( self.num_activesite_residues )
                                                  , 3 ) )
        self.percentage_activesite_HIS = ( round( float( self.AS_HIS ) / 
                                                  float( self.num_activesite_residues )
                                                  , 3 ) )
        self.percentage_activesite_ILE = ( round( float( self.AS_ILE ) / 
                                                  float( self.num_activesite_residues )
                                                  , 3 ) )
        self.percentage_activesite_LYS = ( round( float( self.AS_LYS ) / 
                                                  float( self.num_activesite_residues )
                                                  , 3 ) )
        self.percentage_activesite_LEU = ( round( float( self.AS_LEU ) / 
                                                  float( self.num_activesite_residues )
                                                  , 3 ) )
        self.percentage_activesite_MET = ( round( float( self.AS_MET ) / 
                                                  float( self.num_activesite_residues )
                                                  , 3 ) )
        self.percentage_activesite_ASN = ( round( float( self.AS_ASN ) / 
                                                  float( self.num_activesite_residues )
                                                  , 3 ) )
        self.percentage_activesite_PRO = ( round( float( self.AS_PRO ) / 
                                                  float( self.num_activesite_residues )
                                                  , 3 ) )
        self.percentage_activesite_GLN = ( round( float( self.AS_GLN ) / 
                                                  float( self.num_activesite_residues )
                                                  , 3 ) )
        self.percentage_activesite_ARG = ( round( float( self.AS_ARG ) / 
                                                  float( self.num_activesite_residues )
                                                  , 3 ) )
        self.percentage_activesite_SER = ( round( float( self.AS_SER ) / 
                                                  float( self.num_activesite_residues )
                                                  , 3 ) )
        self.percentage_activesite_THR = ( round( float( self.AS_THR ) / 
                                                  float( self.num_activesite_residues )
                                                  , 3 ) )
        self.percentage_activesite_VAL = ( round( float( self.AS_VAL ) / 
                                                  float( self.num_activesite_residues )
                                                  , 3 ) )
        self.percentage_activesite_TRP = ( round( float( self.AS_TRP ) / 
                                                  float( self.num_activesite_residues )
                                                  , 3 ) )
        self.percentage_activesite_TYR = ( round( float( self.AS_TYR ) / 
                                                  float( self.num_activesite_residues )
                                                  , 3 ) )


        return True



    def get_activesite_AA_composition_per_lig_res( self ):
        self.AS_activesite_res_per_lig_list = []
        
        # for each ligand
        for uniq_lig_name in self.activesite_residues_per_lig.keys():
            # get the unique protein names of the residues around this ligand
            activesite_uniq_pro_res_names = self.activesite_residues_per_lig[ uniq_lig_name ]
            
            # get the total number of amino acids in this activesite
            tot_AA_in_activesite = len( activesite_uniq_pro_res_names )
            
            # now just collect the three-letter amino acid codes of the activesite residues
            activesite_pro_residues = []
            for uniq_pro_name in activesite_uniq_pro_res_names:
                res_name = unqi_pro_name.split( '_' )[0]
                activesite_pro_residues.append( res_name )
            
            # count the number of amino acid residues around the ligand and append to data lists
            ALA_per_lig = ( activesite_pro_residues.count( "ALA" ) )
            CYS_per_lig = ( activesite_pro_residues.count( "CYS" ) )
            ASP_per_lig = ( activesite_pro_residues.count( "ASP" ) )
            GLU_per_lig = ( activesite_pro_residues.count( "GLU" ) )
            PHE_per_lig = ( activesite_pro_residues.count( "PHE" ) )
            GLY_per_lig = ( activesite_pro_residues.count( "GLY" ) )
            HIS_per_lig = ( activesite_pro_residues.count( "HIS" ) )
            ILE_per_lig = ( activesite_pro_residues.count( "ILE" ) )
            LYS_per_lig = ( activesite_pro_residues.count( "LYS" ) )
            LEU_per_lig = ( activesite_pro_residues.count( "LEU" ) )
            MET_per_lig = ( activesite_pro_residues.count( "MET" ) )
            ASN_per_lig = ( activesite_pro_residues.count( "ASN" ) )
            PRO_per_lig = ( activesite_pro_residues.count( "PRO" ) )
            GLN_per_lig = ( activesite_pro_residues.count( "GLN" ) )
            ARG_per_lig = ( activesite_pro_residues.count( "ARG" ) )
            SER_per_lig = ( activesite_pro_residues.count( "SER" ) )
            THR_per_lig = ( activesite_pro_residues.count( "THR" ) )
            VAL_per_lig = ( activesite_pro_residues.count( "VAL" ) )
            TRP_per_lig = ( activesite_pro_residues.count( "TRP" ) )
            TYR_per_lig = ( activesite_pro_residues.count( "TYR" ) )
            
            # collect percentage data of number of a specific amino acid in the 
            # ligand's activesite versus the number of total amino acids in that activesite
            
            
            
            # recollect the percentage of each residue in the protein
            '''
            percentage_tot_per_lig_ALA = ( round( float( ALA_per_lig ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
            percentage_tot_per_lig_CYS = ( round( float( CYS_per_lig ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
            percentage_tot_per_lig_ASP = ( round( float( ASP_per_lig ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
            percentage_tot_per_lig_GLU = ( round( float( GLU_per_lig ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
            percentage_tot_per_lig_PHE = ( round( float( PHE_per_lig ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
            percentage_tot_per_lig_GLY = ( round( float( GLY_per_lig ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
            percentage_tot_per_lig_HIS = ( round( float( HIS_per_lig ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
            percentage_tot_per_lig_ILE = ( round( float( ILE_per_lig ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
            percentage_tot_per_lig_LYS = ( round( float( LYS_per_lig ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
            percentage_tot_per_lig_LEU = ( round( float( LEU_per_lig ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
            percentage_tot_per_lig_MET = ( round( float( MET_per_lig ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
            percentage_tot_per_lig_ASN = ( round( float( ASN_per_lig ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
            percentage_tot_per_lig_PRO = ( round( float( PRO_per_lig ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
            percentage_tot_per_lig_GLN = ( round( float( GLN_per_lig ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
            percentage_tot_per_lig_ARG = ( round( float( ARG_per_lig ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
            percentage_tot_per_lig_SER = ( round( float( SER_per_lig ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
            percentage_tot_per_lig_THR = ( round( float( THR_per_lig ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
            percentage_tot_per_lig_VAL = ( round( float( VAL_per_lig ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
            percentage_tot_per_lig_TRP = ( round( float( TRP_per_lig ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
            percentage_tot_per_lig_TYR = ( round( float( TYR_per_lig ) / 
                                           float( self.size_of_protein )
                                           , 3 ) )
            '''
                                    
        return True



    def collect_activesite_data( self ):
        # append all of the final data to the self.lists
        # because if the analysis got this far, that means there actually is data to collect for the PDB
        self.AS_pdb_names.append( self.name )
        self.AS_lig_res.append( self.num_lig_residues )
        self.AS_lig_atoms.append( self.num_lig_atoms )
        self.AS_activesite_res.append( self.num_activesite_res )
        self.AS_activesite_atoms.append( self.num_activesite_atoms )
        
        self.AS_ALA_list.append( self.AS_ALA )
        self.AS_CYS_list.append( self.AS_CYS )
        self.AS_ASP_list.append( self.AS_ASP )
        self.AS_GLU_list.append( self.AS_GLU )
        self.AS_PHE_list.append( self.AS_PHE )
        self.AS_GLY_list.append( self.AS_GLY )
        self.AS_HIS_list.append( self.AS_HIS )
        self.AS_ILE_list.append( self.AS_ILE )
        self.AS_LYS_list.append( self.AS_LYS )
        self.AS_LEU_list.append( self.AS_LEU )
        self.AS_MET_list.append( self.AS_MET )
        self.AS_ASN_list.append( self.AS_ASN )
        self.AS_PRO_list.append( self.AS_PRO )
        self.AS_GLN_list.append( self.AS_GLN )
        self.AS_ARG_list.append( self.AS_ARG )
        self.AS_SER_list.append( self.AS_SER )
        self.AS_THR_list.append( self.AS_THR )
        self.AS_VAL_list.append( self.AS_VAL )
        self.AS_TRP_list.append( self.AS_TRP )
        self.AS_TYR_list.append( self.AS_TYR )
        
        self.tot_ALA_list.append( self.tot_ALA )
        self.tot_CYS_list.append( self.tot_CYS )
        self.tot_ASP_list.append( self.tot_ASP )
        self.tot_GLU_list.append( self.tot_GLU )
        self.tot_PHE_list.append( self.tot_PHE )
        self.tot_GLY_list.append( self.tot_GLY )
        self.tot_HIS_list.append( self.tot_HIS )
        self.tot_ILE_list.append( self.tot_ILE )
        self.tot_LYS_list.append( self.tot_LYS )
        self.tot_LEU_list.append( self.tot_LEU )
        self.tot_MET_list.append( self.tot_MET )
        self.tot_ASN_list.append( self.tot_ASN )
        self.tot_PRO_list.append( self.tot_PRO )
        self.tot_GLN_list.append( self.tot_GLN )
        self.tot_ARG_list.append( self.tot_ARG )
        self.tot_SER_list.append( self.tot_SER )
        self.tot_THR_list.append( self.tot_THR )
        self.tot_VAL_list.append( self.tot_VAL )
        self.tot_TRP_list.append( self.tot_TRP )
        self.tot_TYR_list.append( self.tot_TYR )
        
        self.AS_num_activesite_nonpolar_atom_types_list.append( self.num_activesite_nonpolar_atoms )
        self.AS_num_activesite_polar_atom_types_list.append( self.num_activesite_polar_atoms )
        self.AS_num_activesite_unk_atom_types_list.append( self.num_activesite_unk_atom_types )
        self.AS_num_ligand_nonpolar_atom_types_list.append( self.num_lig_nonpolar_atoms )
        self.AS_num_ligand_polar_atom_types_list.append( self.num_lig_polar_atoms )
        self.AS_num_ligand_metal_atom_types_list.append( self.num_lig_metal_atoms )
        self.AS_num_ligand_unk_atom_types_list.append( self.num_lig_unk_atoms )
        
        self.percentage_activesite_nonpolar_list.append( self.percentage_activesite_nonpolar )
        self.percentage_activesite_polar_list.append( self.percentage_activesite_polar )
        self.percentage_ligand_nonpolar_list.append( self.percentage_ligand_nonpolar )
        self.percentage_ligand_polar_list.append( self.percentage_ligand_polar )
        self.percentage_ligand_metal_list.append( self.percentage_ligand_metal )
        
        self.percentage_activesite_ALA_list.append( self.percentage_activesite_ALA )
        self.percentage_activesite_CYS_list.append( self.percentage_activesite_CYS )
        self.percentage_activesite_ASP_list.append( self.percentage_activesite_ASP )
        self.percentage_activesite_GLU_list.append( self.percentage_activesite_GLU )
        self.percentage_activesite_PHE_list.append( self.percentage_activesite_PHE )
        self.percentage_activesite_GLY_list.append( self.percentage_activesite_GLY )
        self.percentage_activesite_HIS_list.append( self.percentage_activesite_HIS )
        self.percentage_activesite_ILE_list.append( self.percentage_activesite_ILE )
        self.percentage_activesite_LYS_list.append( self.percentage_activesite_LYS )
        self.percentage_activesite_LEU_list.append( self.percentage_activesite_LEU )
        self.percentage_activesite_MET_list.append( self.percentage_activesite_MET )
        self.percentage_activesite_ASN_list.append( self.percentage_activesite_ASN )
        self.percentage_activesite_PRO_list.append( self.percentage_activesite_PRO )
        self.percentage_activesite_GLN_list.append( self.percentage_activesite_GLN )
        self.percentage_activesite_ARG_list.append( self.percentage_activesite_ARG )
        self.percentage_activesite_SER_list.append( self.percentage_activesite_SER )
        self.percentage_activesite_THR_list.append( self.percentage_activesite_THR )
        self.percentage_activesite_VAL_list.append( self.percentage_activesite_VAL )
        self.percentage_activesite_TRP_list.append( self.percentage_activesite_TYP )
        self.percentage_activesite_TYR_list.append( self.percentage_activesite_TYR )
        
        self.percentage_tot_ALA_list.append( self.percentage_tot_ALA )
        self.percentage_tot_CYS_list.append( self.percentage_tot_CYS )
        self.percentage_tot_ASP_list.append( self.percentage_tot_ASP )
        self.percentage_tot_GLU_list.append( self.percentage_tot_GLU )
        self.percentage_tot_PHE_list.append( self.percentage_tot_PHE )
        self.percentage_tot_GLY_list.append( self.percentage_tot_GLY )
        self.percentage_tot_HIS_list.append( self.percentage_tot_HIS )
        self.percentage_tot_ILE_list.append( self.percentage_tot_ILE )
        self.percentage_tot_LYS_list.append( self.percentage_tot_LYS )
        self.percentage_tot_LEU_list.append( self.percentage_tot_LEU )
        self.percentage_tot_MET_list.append( self.percentage_tot_MET )
        self.percentage_tot_ASN_list.append( self.percentage_tot_ASN )
        self.percentage_tot_PRO_list.append( self.percentage_tot_PRO )
        self.percentage_tot_GLN_list.append( self.percentage_tot_GLN )
        self.percentage_tot_ARG_list.append( self.percentage_tot_ARG )
        self.percentage_tot_SER_list.append( self.percentage_tot_SER )
        self.percentage_tot_THR_list.append( self.percentage_tot_THR )
        self.percentage_tot_VAL_list.append( self.percentage_tot_VAL )
        self.percentage_tot_TRP_list.append( self.percentage_tot_TRP )
        self.percentage_tot_TYR_list.append( self.percentage_tot_TYR )
        
        return True
