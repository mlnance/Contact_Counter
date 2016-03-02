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


###############################
#### CONTACT COUNTING CODE ####
###############################

class CTCT:
    def __init__( self ):
        # get current working directory
        self.working_dir = os.getcwd() + '/'
        
        # make data lists to add over course of program for contact counts
        # will be added to pandas df at end
        self.CC_pdb_names = []
        self.CC_lig_atms = []
        self.CC_activesite_atms = []
        self.CC_pp_contacts = []
        self.CC_pn_contacts = []
        self.CC_np_contacts = []
        self.CC_nn_contacts = []
        self.CC_unk_contacts = []
        
        # make data lists to add over course of program for contact counts per lig
        # will be added to pandas df at end
        self.CC_per_lig_pdb_names = []
        self.CC_per_lig_lig_names = []
        self.CC_per_lig_lig_atms = []
        self.CC_per_lig_activesite_atms = []
        self.CC_per_lig_pp_contacts = []
        self.CC_per_lig_pn_contacts = []
        self.CC_per_lig_np_contacts = []
        self.CC_per_lig_nn_contacts = []
        self.CC_per_lig_unk_contacts = []
        
        

    def read_pro_lig_pickles( self, pro_pickle, lig_pickle ):
        # get the four letter code from the passed protein pickle name
        self.name = pro_pickle.split( '/' )[-1][:4].lower()
        
        # collect the corresponding protein and ligand dictionary data
        self.protein = pickle.load( open( pro_pickle, "rb" ) )
        self.ligand = pickle.load( open( lig_pickle, "rb" ) )
        
        # collect number of ligand residues and atoms
        self.num_lig_res = len( self.ligand.keys() )
        self.num_lig_atoms = 0
        for lig_res in self.ligand.keys():
            self.num_lig_atoms += len( self.ligand[ lig_res ] )
        
        return True
    
    
    
    def get_activesite( self, cutoff ):
        # overall activesite dictionary
        # key: unique protein name (resname_reschain_resnum), value: list of ATOM lines per residue
        self.activesite_dict = {}
        
        # unique ligand name (resname_reschain_resnum)
        #value: list of 3-letter amino acid names
        self.activesite_lig_pro_res_dict = {}
        
        # unique ligand name (resname_reschain_resnum)
        #value: list of atom_line for each AA ( to get atom count later )
        self.activesite_lig_pro_atoms_dict = {}
        
        # activesite data holders
        self.num_activesite_res = 0
        self.activesite_residues = []
        self.num_activesite_atms = 0
        self.num_activesite_nonpolar_atoms = 0
        self.num_activesite_polar_atoms = 0
        self.num_activesite_unk_atom_types = 0
        self.activesite_num_nonpolar_atoms = {}
        self.activesite_num_polar_atoms = {}
        self.activesite_num_unk_atoms = {}
        
        for uniq_lig_name in self.ligand.keys():
            # list to store the 3 letter names of all of the protein residues within the cutoff distance of each ligand residue (by unique name)
            AS_names_in_activesite = []
            AS_atms_in_activesite = []
            
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
                                three_letter_name = uniq_pro_name[0:3]
                                AS_names_in_activesite.append( three_letter_name )
                                
                            # store atom_line for each unique amino acid within the activesite to get an atom count later
                            if atom_line not in AS_atms_in_activesite:
                                AS_atms_in_activesite.append( atom_line )
                        
            # store the list of the 3 letter names for the amino acid within the self.activesite_lig_pro_dict according to which ligand it is near
            self.activesite_lig_pro_res_dict[ uniq_lig_name ] = AS_names_in_activesite
            self.activesite_lig_pro_atoms_dict[ uniq_lig_name ] = AS_atms_in_activesite
            
        # get number of activesite residues
        self.num_activesite_res = len( self.activesite_dict.keys() )
        
        # count the number of activesite atoms
        for uniq_lig_name in self.activesite_lig_pro_atoms_dict.keys():
            self.num_activesite_atms += len( self.activesite_lig_pro_atoms_dict[ uniq_lig_name ] )
            
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
            print "and", self.num_activesite_atms, "non-hydrogen activesite atoms"
            return True
        

        
    def count_contacts( self, cutoff ):
        # must have already found all residues within the activesite
        # counted as ligand to protein!!
        self.polar_polar = 0
        self.polar_nonpolar = 0
        self.nonpolar_polar = 0
        self.nonpolar_nonpolar = 0
        self.unk_contact = 0
        
        # ligxyz_proxyz xyz coordinates unique for every atom, best way to collect unique contacts made
        self.uniq_contact_list = []
        
        for uniq_lig_name in self.ligand.keys():
            # add pdb name and ligand name to list
            self.CC_per_lig_pdb_names.append( self.name )
            self.CC_per_lig_lig_names.append( uniq_lig_name )
            self.CC_per_lig_lig_atms.append( len( self.ligand[ uniq_lig_name ] ) )
            self.CC_per_lig_activesite_atms.append( len( self.activesite_lig_pro_atoms_dict[ uniq_lig_name ] ) )
            
            # make empty counters - used for counting contacts per ligand
            lig_polar_polar = 0
            lig_polar_nonpolar = 0
            lig_nonpolar_polar = 0
            lig_nonpolar_nonpolar = 0
            lig_unk_contact = 0
            
            # for every ligand residue
            for lig_pdb_line in self.ligand[ uniq_lig_name ]:
                # get ligand atom xyz
                x_lig = lig_pdb_line.x_coord
                y_lig = lig_pdb_line.y_coord
                z_lig = lig_pdb_line.z_coord
                lig_xyz = [ x_lig, y_lig, z_lig ]
                lig_xyz_str = str( x_lig ) + '_' + str( y_lig ) + '_' + str( z_lig )
                
                for pro_pdb_line in self.activesite_lig_pro_atoms_dict[ uniq_lig_name ]:
                    # get protein atom xyz
                    x_pro = pro_pdb_line.x_coord
                    y_pro = pro_pdb_line.y_coord
                    z_pro = pro_pdb_line.z_coord
                    pro_xyz = [ x_pro, y_pro, z_pro ]
                    pro_xyz_str = str( x_pro ) + '_' + str( y_pro ) + '_' + str( z_pro )
                    
                    # check atomic distance
                    contact_distance = calc_distance( lig_xyz, pro_xyz )
                    if contact_distance < cutoff:
                        # check to see that this contact has not yet been counted
                        uniq_contact = lig_xyz_str + '.' + pro_xyz_str
                        if uniq_contact not in self.uniq_contact_list:
                            # append unique contact to list
                            self.uniq_contact_list.append( uniq_contact )
                            
                            # check contact type
                            if lig_pdb_line.element in polar_atoms or lig_pdb_line.element in metal_atoms:
                                # polar polar
                                if pro_pdb_line.element in polar_atoms or pro_pdb_line.element in metal_atoms:
                                    lig_polar_polar += 1
                                    
                                # polar nonpolar
                                elif pro_pdb_line.element in nonpolar_atoms:
                                    lig_polar_nonpolar += 1
                                   
                                # unknown
                                else:
                                    lig_unk_contact += 1
                                   
                            elif lig_pdb_line.element in nonpolar_atoms:
                                # nonpolar polar
                                if pro_pdb_line.element in polar_atoms or pro_pdb_line.element in metal_atoms:
                                    lig_nonpolar_polar += 1
                                   
                                # nonpolar nonpolar
                                elif pro_pdb_line.element in nonpolar_atoms:
                                    lig_nonpolar_nonpolar += 1
                                    
                                # unknown
                                else:
                                    lig_unk_contact += 1
                                   
                            # both unknown
                            else:
                                lig_unk_contact += 1
                                
            # append this info to the per ligand data frame lists
            self.CC_per_lig_pp_contacts.append( lig_polar_polar )
            self.CC_per_lig_pn_contacts.append( lig_polar_nonpolar )
            self.CC_per_lig_np_contacts.append( lig_nonpolar_polar )
            self.CC_per_lig_nn_contacts.append( lig_nonpolar_nonpolar )
            self.CC_per_lig_unk_contacts.append( lig_unk_contact )
            
            # add the counts for this particular ligand to the total count
            self.polar_polar += lig_polar_polar
            self.polar_nonpolar += lig_polar_nonpolar
            self.nonpolar_polar += lig_nonpolar_polar
            self.nonpolar_nonpolar += lig_nonpolar_nonpolar
            self.unk_contact += lig_unk_contact
            
            
        # store all data in global list
        # because if it got here, all of the data existed
        self.CC_pdb_names.append( self.name )
        self.CC_lig_atms.append( self.num_lig_atoms )
        self.CC_activesite_atms.append( self.num_activesite_atms )
        self.CC_pp_contacts.append( self.polar_polar )
        self.CC_pn_contacts.append( self.polar_nonpolar )
        self.CC_np_contacts.append( self.nonpolar_polar )
        self.CC_nn_contacts.append( self.nonpolar_nonpolar )
        self.CC_unk_contacts.append( self.unk_contact )
        
        return True
