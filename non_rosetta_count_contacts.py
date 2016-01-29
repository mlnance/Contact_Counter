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
try:
    from colorama import Fore, Style
except:
    pass
try:
    from Bio.PDB import *
except ImportError:
    print "Trouble with imports - do you have Bio.PDB? Exiting"
    sys.exit()
try:
    import networkx as nx
except ImportError:
    print "Trouble with imports - do you have networkx? Exiting"
    sys.exit()



#######################
#### CHEMICAL DATA ####
#######################

# list of polar and nonpolar atoms
nonpolar_atoms = [ 'C' ]
polar_atoms = [ 'O', 'N', 'S', 'P', 'F', 'Cl', 'CL', 'Br', 'BR', 'I', 'Se', 'SE', 'B' ]

# list of residues to remove if they are a ligand
# ligand meaning it was defined as HETATM in the PDB
metal_list = [ 'B', 'K', 'V', 'Y', 'W', 'U', "LI", "BE", "NA", "MG", "AL", "SI", "CA", "SC", "TI", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "GE", "AS", "RB", "SR", "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN", "SB", "TE", "CS", "BA", "LU", "HF", "TA", "RE", "OS", "IR", "PT", "AU", "HG", "TI", "PB", "BI", "PO", "AT", "FR", "RA", "LR", "RF", "DB", "SG", "BH", "HS", "MT", "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "AC", "TH", "PA", "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM", "MD", "NO" ]

AA_list = [ "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TRY", "VAL" ]

nucleic_acid_list = [ 'DA', 'A', 'DT', 'T', 'DG', 'G', 'DC', 'C' ]



##############################
#### PDB LINE DEFINITIONS ####
##############################

class PDB_line:
    def __init__(self, line):
        # only works for ATOM and HETATM lines ( TER lines don't have coordinates for instance, thus would crash )
        # line locations taken from the PDB databases's pdb format description page
        # ftp://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf
        # line 187 for ATOM
        # line 194 for HETATM
        self.line = line.rstrip( '/n' )
        
    def atom_num(self):
        """
        Returns the atom number from a PDB ATOM line.
        Removes any leftover white space.
        Example: '1', "1234"
        :return: int( atom number )
        """
        return int( self.line[6:11].replace( ' ', '' ) )
    
    def atom_name(self):
        """
        Returns the atom name from a PDB ATOM line.
        Removes any leftover white space.
        Example: 'C', "Br"
        :return: str( atom name )
        """
        return str( self.line[12:16].replace( ' ', '' ) )
    
    def alt_loc(self):
        """
        Returns the alternate location identifier from a PDB ATOM line.
        Example: 'A' or 'B'
        :return: str( alternate location identifier )
        """
        return str( self.line[16:17].replace( ' ', '' ) )
    
    def res_name(self):
        """
        Returns the residue name from a PDB ATOM line.
        Removes any leftover white space.
        Example: "GLY", "NAG"
        :return: str( residue name )
        """
        return str( self.line[17:20].replace( ' ', '' ) )
    
    def res_chain(self):
        """
        Returns the residue chain identifier from a PDB ATOM line.
        Removes any leftover white space.
        Example: 'A', 'B', 'C'
        :return: str( residue name )
        """
        return str( self.line[21:22] )
    
    def res_num(self):
        """
        Returns the residue number from a PDB ATOM line.
        Removes any leftover white space.
        Example: '4', "1000"
        :return: int( residue number )
        """
        return int( self.line[22:26].replace( ' ', '' ) )
    
    def i_code(self):
        """
        Returns the insertion code from a PDB ATOM line.
        :return: str( insertion code )
        """
        return str( self.line[26:27] )
    
    def res_x_coord(self):
        """
        Returns the x-coordinate from a PDB ATOM line.
        :return: float( x-cooridnate)
        """
        return float( self.line[30:38].replace( ' ', '' ) )
    
    def res_y_coord(self):
        """
        Returns the y-coordinate from a PDB ATOM line.
        :return: float( y-cooridnate)
        """
        return float( self.line[38:46].replace( ' ', '' ) )
    
    def res_z_coord(self):
        """
        Returns the z-coordinate from a PDB ATOM line.
        :return: float( z-cooridnate)
        """
        return float( self.line[46:54].replace( ' ', '' ) )
    
    def occupancy(self):
        """
        Returns the occupancy value from a PDB ATOM line.
        :return: float( occupancy )
        """
        return float( self.line[54:60].replace( ' ', '' ) )

    def temp_factor(self):
        """
        Returns the temperature factor from a PDB ATOM line.
        :return: float( temperature factor )
        """
        return float( self.line[60:66].replace( ' ', '' ) )
    
    def element(self):
        """
        Returns the one- or two-letter atomic element name from a PDB ATOM line.
        :return: str( element )
        """
        return str( self.line[76:78].replace( ' ', '' ) )
    
    def charge(self):
        """
        Returns the atomic charge from a PDB ATOM line.
        :return: str( atomic charge )
        """
        return str( self.line[78:80].replace( ' ', '' ) )



class LINK_line:
    def __init__(self, line):
        # only works for LINK lines
        # line locations taken from the PDB databases's pdb format description page
        # ftp://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf
        # page 173
        self.line = line.rstrip( '/n' )
        
    def atom_name1(self):
        return str( self.line[ 12:16 ] )
    
    def alt_loc1(self):
        return str( self.line[ 16:17] )
    
    def res1_name(self):
        return str( self.line[17:20] )
    
    def res1_chain(self):
        return str( self.line[21:22] )
    
    def res1_seq(self):
        return int( self.line[22:26] )
    
    def i_code1(self):
        return str( self.line[26:27] )
    
    def atom_name2(self):
        return str( self.line[42:46] )
    
    def alt_loc2(self):
        return str( self.line[46:47] )
    
    def res2_name(self):
        return str( self.line[47:50] )
    
    def res2_chain(self):
        return str( self.line[51:52] )
    
    def res2_seq(self):
        return int( self.line[52:56] )
    
    def i_code2(self):
        return str( self.line[56:57] )
    
    def link_dist(self):
        return int( self.line[73:78] )



class MODRES_line:
    def __init__(self, line):
        # only works for MODRES lines
        # line locations taken from the PDB databases's pdb format description page
        # ftp://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf
        # page 157
        self.line = line.rstrip( '/n' )

    def id_code(self):
        return int( self.line[7:11].replace( ' ', '' ) )

    def res_name(self):
        # the residue name used in the PDB file (post-modification)
        return str( self.line[12:15].replace( ' ', '' ) )
    
    def res_chain(self):
        return str( self.line[16:17] )
    
    def res_num(self):
        return int( self.line[18:22].replace( ' ', '' ) )
    
    def i_code(self):
        # insertion code
        return str( self.line[22:23] )
    
    def std_res_name(self):
        # what the standard residue name is (pre-modification)
        return str( self.line[24:27].replace( ' ', '' ) )
    
    def comment(self):
        return str( self.line[29:70] )



##############################
#### CIF LINE DEFINITIONS ####
##############################

class cif_struct_conn_lines:
    def __init__(self, cif_filename):
        """
        Give this class an object (list) of the lines in the cif file
        """
        # only works for the struct_conn line in an mmcif file (ie a xxxx.cif file)
        self.cif_filename = cif_filename
    
    def check_if_has_mmcif_dict(self):
        try:
            self.mmcif_dict = MMCIF2Dict.MMCIF2Dict( self.cif_filename )
            self.mmcif_dict[ "_struct_conn.conn_type_id" ]
            self.mmcif_dict[ "_struct_conn.ptnr1_auth_asym_id" ] 
            self.mmcif_dict[ "_struct_conn.ptnr1_auth_comp_id" ]
            self.mmcif_dict[ "_struct_conn.ptnr1_auth_seq_id" ]
            self.mmcif_dict[ "_struct_conn.ptnr2_auth_asym_id" ] 
            self.mmcif_dict[ "_struct_conn.ptnr2_auth_comp_id" ]
            self.mmcif_dict[ "_struct_conn.ptnr2_auth_seq_id" ]
            return True
        except:
            # print in color (if possible) letting the user know the cif file didn't work
            try:
                text = "~~Something wrong with the _struct_conn keys in " + self.cif_filename.split( '/' )[-1] + " using LINK records instead"
                print(Fore.YELLOW + text + Style.RESET_ALL)
            except:
                print "~~Something wrong with the _struct_conn keys in", self.cif_filename.split( '/' )[-1], "using LINK records instead"
            return False
        
    def connection_types(self):
        return self.mmcif_dict[ "_struct_conn.conn_type_id" ]
    
    def partner1_chain_ids(self):
        return self.mmcif_dict[ "_struct_conn.ptnr1_auth_asym_id" ] 
   
    def partner1_res_names(self):
        return self.mmcif_dict[ "_struct_conn.ptnr1_auth_comp_id" ]
    
    def partner1_res_nums(self):
        return self.mmcif_dict[ "_struct_conn.ptnr1_auth_seq_id" ]
    
    def partner2_chain_ids(self):
        return self.mmcif_dict[ "_struct_conn.ptnr2_auth_asym_id" ]
    
    def partner2_res_names(self):
        return self.mmcif_dict[ "_struct_conn.ptnr2_auth_comp_id" ]
    
    def partner2_res_nums(self):
        return self.mmcif_dict[ "_struct_conn.ptnr2_auth_seq_id" ]
    
    def get_uniq_connection_names(self):
        connection_types = self.connection_types()
        res1_names = self.partner1_res_names()
        res1_chains = self.partner1_chain_ids()
        res1_nums = self.partner1_res_nums()
        res2_names = self.partner2_res_names()
        res2_chains = self.partner2_chain_ids()
        res2_nums = self.partner2_res_nums()
        
        unique_connection_names = []
        for ii in range( len( connection_types ) ):
            if connection_types[ii] == "covale":
                res1_unique_name = res1_names[ii] + '_' + res1_chains[ii] + '_' + res1_nums[ii]
                res2_unique_name = res2_names[ii] + '_' + res2_chains[ii] + '_' + res2_nums[ii]
                uniq_connection_name = res1_unique_name + '+' + res2_unique_name
                if uniq_connection_name not in unique_connection_names:
                    unique_connection_names.append( uniq_connection_name )
        return unique_connection_names        
    


class cif_atom_site_lines:
    def __init__(self, cif_filename):
        """
        Give this class an object (list) of the lines in the cif file
        """
        # only works for the atom_site line in an mmcif file (ie a xxx.cif file)
        try:
            self.cif_filename = cif_filename
            self.mmcif_dict = MMCIF2Dict.MMCIF2Dict( self.cif_filename )
            self.mmcif_dict[ "_atom_site.group_PDB" ]
            self.mmcif_dict[ "_atom_site.auth_asym_id" ]
            self.mmcif_dict[ "_atom_site.auth_comp_id" ]
            self.mmcif_dict[ "_atom_site.auth_seq_id" ]
        except:
            pass
        
    def atom_types(self):
        """
        Returns a list of "ATOM" and "HETATM" strings
        """
        return self.mmcif_dict[ "_atom_site.group_PDB" ]
    
    def chain_ids(self):
        return self.mmcif_dict[ "_atom_site.auth_asym_id" ]
    
    def res_names(self):
        return self.mmcif_dict[ "_atom_site.auth_comp_id" ]
    
    def res_nums(self):
        return self.mmcif_dict[ "_atom_site.auth_seq_id" ]
    
    def cif_uniq_pro_names(self):
        atom_types = self.atom_types()
        res_names = self.res_names()
        res_chains = self.chain_ids()
        res_nums = self.res_nums()
        
        unique_protein_names = []
        for ii in range( len( atom_types ) ):
            if atom_types[ii] == "ATOM":
                unique_name = res_names[ii] + '_' + res_chains[ii] + '_' + res_nums[ii]
                if unique_name not in unique_protein_names:
                    unique_protein_names.append( unique_name )
        return unique_protein_names

    def cif_uniq_het_names(self):
        atom_types = self.atom_types()
        res_names = self.res_names()
        res_chains = self.chain_ids()
        res_nums = self.res_nums()
        
        unique_hetatm_names = []
        for ii in range( len( atom_types ) ):
            if atom_types[ii] == "HETATM":
                if res_names[ii] != "HOH" and res_names[ii] != "DOD":
                    unique_name = res_names[ii] + '_' + res_chains[ii] + '_' + res_nums[ii]
                    if unique_name not in unique_hetatm_names:
                        unique_hetatm_names.append( unique_name )
        return unique_hetatm_names



###############################
#### CONTACT COUNTING CODE ####
###############################

class CTCT:
    def __init__( self, pdb_name_list, download_pdbs ):
        """
        Initializes variables that hold the data, and reads the PDB name list read in on startup
        """
        # get current working directory
        self.working_dir = os.getcwd()
        
        # load up PDB name list from input argument
        self.pdb_names = []
        pdb_file_CR = open( pdb_name_list, 'r' ).readlines()  # CR = carriage return (new line character)
        for name in pdb_file_CR:
            name = name.rstrip( '\n' )
            # add .pdb extension if not already there  -  don't need to do this if downloading the PDB
            if not download_pdbs:
                if not name.endswith( ".pdb" ):
                    name = name + ".pdb"
            self.pdb_names.append( name )

    
    def instantiate_holders( self ):
        # instantiate lists that will hold relevant PDB data
        self.protein_lines = []
        self.hetatm_lines = []
        self.link_records = []
        self.protein = {}
        self.ligand = {}


    def instantiate_data_holders( self ):
        ## initialize all the data holders
        # make on/off switches for skipping 
        self.glycosylated_proteins = []
        
        # hold the name of PDBs that contain undesirables like multiple models or an unknown HETATM
        self.unknown_res_pdb_names = []
        self.multiple_models_pdb_names = []
        self.deuterium_pdb_names = []
        
        # make data lists to add over course of program for AS composition
        self.AS_pdb_names = []
        self.AS_lig_res = []
        self.AS_lig_atms = []
        self.AS_activesite_res = []
        self.AS_activesite_atms = []
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


        # make data lists to add over course of program for AA composition per ligand residue in each pdb
        # each pdb name should show up as many times as it has ligand residues that fit the user's criteria
        self.AS_pdb_names_per_lig = []
        self.AS_lig_uniq_res_names_per_lig = []
        self.AS_lig_res_names_per_lig = []
        self.AS_lig_atms_per_lig = []
        self.AS_activesite_res_per_lig = []
        self.AS_activesite_atms_per_lig = []
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


        # make data lists to add over course of program for contact counts - will be added to pandas df at end
        self.CC_pdb_names = []
        self.CC_lig_atms = []
        self.CC_activesite_atms = []
        self.CC_pp_contacts = []
        self.CC_pn_contacts = []
        self.CC_np_contacts = []
        self.CC_nn_contacts = []
        self.CC_unk_contacts = []


        # make data lists to add over course of program for contact counts per lig - will be added to pandas df at end
        self.CC_per_lig_pdb_names = []
        self.CC_per_lig_lig_names = []
        self.CC_per_lig_lig_atms = []
        self.CC_per_lig_activesite_atms = []
        self.CC_per_lig_pp_contacts = []
        self.CC_per_lig_pn_contacts = []
        self.CC_per_lig_np_contacts = []
        self.CC_per_lig_nn_contacts = []
        self.CC_per_lig_unk_contacts = []



    def download_pdb( self, pdb_name ):
        """
        Uses a utility function to download a PDB from the internet based on the four letter accession code
        :param pdb_name: str( four-letter PDB code )
        :return: str( full path name to the PDB file )
        """
        from util import download_pdb
        
        # relevant paths
        cur_dir = os.getcwd() + '/'
        pdb_dir = cur_dir + 'pdbs/'
        
        # check to make sure the 'pdbs' directory exists one down from the current directory
        if not os.path.isdir( pdb_dir ):
            os.mkdir( pdb_dir )
            
        # if the PDB exists in the 'pdbs' directory, return the path to it
        if os.path.isfile( pdb_dir + pdb_name ):
            return pdb_dir + pdb_name
        
        # otherwise download the PDB
        else:
            # download PDB and get its filename
            try:
                pdb = download_pdb( pdb_name[:4], pdb_dir )
            except:
                pass

            # return the full path name to the PDB file
            return pdb_dir + pdb



    def download_cif( self, pdb_name ):
        """
        Uses a utility function to download the cif file of a PDB from the internet based on the four letter accession code
        :param pdb_name: str( four-letter PDB code )
        :return: str( full path name to the .cif file )
        """
        from util import download_cif_file
        
        # relevant paths
        cur_dir = os.getcwd() + '/'
        pdb_dir = cur_dir + 'pdbs/'
        
        # check to make sure the 'pdbs' directory exists one down from the current directory
        if not os.path.isdir( pdb_dir ):
            os.mkdir( pdb_dir )
            
        # if the PDB exists in the 'pdbs' directory, return the path to it
        if os.path.isfile( pdb_dir + pdb_name ):
            return pdb_dir + pdb_name
        
        # otherwise download the PDB
        else:
            # download PDB and get its filename
            try:
                pdb = download_cif_file( pdb_name[:4], pdb_dir )
            except:
                pass
            
            # return the full path name to the PDB file
            return pdb_dir + pdb



    def pymol_clean( self, pdb_filename ):
        """
        Uses PyMOL to "clean" a PDB by removing waters and hydrogens.
        Opens up PyMOL via a command-line argument. Closes PyMOL upon completion
        :param pdb_filename:
        Doesn't return anything. Just saves the cleaned PDB file
        """
        # get the right names of the PDB files
        split_pdb_name = pdb_filename.split( '/' )[-1]
        pdb_name = split_pdb_name[:-4]
        
        # use pymol to remove water residues and remove hydrogens  -  ie. clean the file
        pymol_command = "%s -cqd 'load %s; remove resn hoh; remove resn dod; remove hydrogens %s; save %s'" %( pymol_exe_dir, pdb_filename, pdb_name, pdb_name )
        
        # run the PyMol command through the terminal
        os.popen( pymol_command )
        
        
        
    def get_uniq_connection_names_from_LINK_records( self, LINK_records ):
        unique_connection_names = []
        
        for link_line in LINK_records:
            res1_unique_name = link_line.res1_name + '_' + link_line.res1_chain + '_' + link_line.res1_seq
            res2_unique_name = link_line.res2_name + '_' + link_line.res2_chain + '_' + link_line.res2_seq
            uniq_connection_name = res1_unique_name + '+' + res2_unique_name
            
            if uniq_connection_name not in unique_connection_names:
                unique_connection_names.append( uniq_connection_name )
        
        return unique_connection_names


    
    def graph_out_residue_connections( self, unique_partner_names, unique_protein_names, unique_hetatm_names ):
        # make a tree using networkx
        graph = nx.Graph()
        
        # add nodes 2 at a time from unique_connect_partner_names
        for partners in unique_partner_names:
            # get the unique names for partner1 and partner2
            ptnr1 = partners.split( '+' )[0]
            ptnr2 = partners.split( '+' )[1]
            
            # get the residue name, chain, and sequence id from the unique names
            ptnr1_resname = ptnr1.split( '_' )[0]
            ptnr1_reschain = ptnr1.split( '_' )[1]
            ptnr1_resseq = ptnr1.split( '_' )[2]
            ptnr2_resname = ptnr2.split( '_' )[0]
            ptnr2_reschain = ptnr2.split( '_' )[1]
            ptnr2_resseq = ptnr2.split( '_' )[2]
            
            ## add each node according to whether its a HETATM or not
            # if residue not already a node
            if not ptnr1 in graph.nodes():
                # add residue as node based on whether it is a protein or hetatm 
                if ptnr1 in unique_protein_names:
                    graph.add_node( ptnr1, { "Chain" : ptnr1_reschain, "Seq" : ptnr1_resseq, "HETATM" : False } )
                else:
                    graph.add_node( ptnr1, { "Chain" : ptnr1_reschain, "Seq" : ptnr1_resseq, "HETATM" : True } )

            # if residue not already a node
            if not ptnr2 in graph.nodes():
                # add residue as node based on whether it is a protein or hetatm
                if ptnr2 in unique_protein_names:
                    graph.add_node( ptnr2, { "Chain" : ptnr2_reschain, "Seq" : ptnr2_resseq, "HETATM" : False } )
                else:
                    graph.add_node( ptnr2, { "Chain" : ptnr2_reschain, "Seq" : ptnr2_resseq, "HETATM" : True } )
                    
            # add a bond (an edge) between the two residues
            graph.add_edge( ptnr1, ptnr2 )
            
        # iter through each subgraph (the group of nodes that have edges connecting them together)
        C = nx.connected_component_subgraphs( graph )
        
        # for holding unique residue names
        remove_these_ligs = []
        
        # for each subgraph
        for g in C:
            # remove = False. Default is we would keep these ligand residues
            remove = False
            for n1, n2 in g.edges_iter():
                # if one node is a HETATM and one is not
                if ( g.node[n1]["HETATM"] == True and g.node[n2]["HETATM"] == False ) or ( g.node[n1]["HETATM"] == False and g.node[n2]["HETATM"] == True ):
                    # remove this subgraph as there is a covalent connection bewteen protein and ligand
                    remove = True
            
                # if you need to remove the ligs from this subgraph
                if remove:
                    for mynode in g.nodes_iter():
                        # if this node is a ligand
                        if g.node[mynode]["HETATM"]:
                            # if you haven't already added it
                            if mynode not in remove_these_ligs:
                                remove_these_ligs.append( mynode )

        # return the unique names of the residues that are covalently linked to the protein
        return remove_these_ligs



    def determine_covalently_bound_ligands( self, pdb_filename, unique_protein_names, unique_hetatm_names, link_records ):
        # download the mmcif file
        cif_filename = self.download_cif( pdb_filename[:4] )
        
        '''
        # get the unique protein and hetatm names from the mmcif file
        # unique name = resname_reschain_resnum
        _atom_site = cif_atom_site_lines( cif_filename )
        unique_protein_names = _atom_site.cif_uniq_pro_names()
        unique_hetatm_names = _atom_site.cif_uniq_het_names()
        '''
        
        # get _struct_conn lines to determine HETATM connections
        _struct_conn = cif_struct_conn_lines( cif_filename )
        
        # check to see if check_if_has_mmcif_dict() by using return value
        response = _struct_conn.check_if_has_mmcif_dict()
        
        # if this PDB has an mmcif file
        if response is True:        
            # unique name = res1name_res1chain_res1num+res2name_res2chain_res2num
            unique_partner_names = _struct_conn.get_uniq_connection_names()
            
            # get list of ligand residues to remove
            remove_these_ligs = self.graph_out_residue_connections( unique_partner_names, unique_protein_names, unique_hetatm_names )
        
        # otherwise this PDB didn't have an mmcif file, use LINK records instead
        else:
            # unique name = res1name_res1chain_res1num+res2name_res2chain_res2num
            unique_partner_names = self.get_uniq_connection_names_from_LINK_records( link_records )
            
            # get list of ligand residues to remove
            remove_these_ligs = self.graph_out_residue_connections( unique_partner_names, unique_protein_names, unique_hetatm_names )

            
        return remove_these_ligs
        
    

    def split_pdb_file( self, pdb_filename, ignore_glycosylated_proteins, keep_clean_pdbs ):
        """
        Splits up the PDB file into ATOM and HETATM lines.
        Doesn't keep HETATM lines that are 1) lone metal atoms, 2) amino acids, or 3) unknown
        Is able to handle modified amino acid residues - potential bug is if the modified residue is a ligand
        Will write out an "error" message if the PDB 1) doesn't have a ligand, 2) has unknown residues, or 3) has multiple models in the file
        :param pdb_filename: str( /path/to/pdb/file )
        :return: 1 if PDB file passes all specifications; 0 otherwise
        """
        # glycosylated protein boolean tracker
        glycosylated_protein = False
        
        # get the PDB name from the end of the full path given in pdb_filename
        split_pdb_name = pdb_filename.split( '/' )[-1]
        pdb_name = split_pdb_name[:-4].lower()
        
        # instantiate lists that will hold relevant PDB data that is NOT wanted
        AA_lig = []
        nuc_acid_lig = []
        water = []
        metals = []
        models = []
        deuterium = []
        unknown = []

        # list of MODRES names to treat as ATOM lines
        modres_protein_res_names = []

        # open up the PDB file
        try:
            with open( pdb_filename, 'r' ) as pdb_fh:
                pdb = pdb_fh.readlines()
        except IOError:
            text = pdb_filename + " doesn't exist in this directory. Did you mean to download it? Exiting."
            print(Fore.BLUE + text + Style.RESET_ALL)
            sys.exit()
        
        # parse through the PDB file
        for line in pdb:
            # if there is a residue that has a modification, it is likely classified as a HETATM, so treat it like an ATOM line
            if line[0:6] == "MODRES":
                modres_line = MODRES_line( line )
                
                # unique residue name
                modres_name = modres_line.res_name()
                modres_chain = modres_line.res_chain()
                modres_num = str( modres_line.res_num() )
                uniq_modres_name = modres_name + '_' + modres_chain + '_' + modres_num
                
                # the standard res name if this residue wasn't modified
                std_res_name = modres_line.std_res_name()
                
                # if the standard residue name is a standard amino acid, add its modified res name to a list to be added later
                if std_res_name in AA_list:
                    modres_protein_res_names.append( uniq_modres_name )
                                
            # if there are multiple models in this PDB file, add the line so this PDB will be skipped
            if line[0:5] == "MODEL":
                models.append( line )
                self.multiple_models_pdb_names.append( pdb_name )
                break
            
            if line[0:4] == "LINK":
                # store LINK records used for determining any covalently bound ligands
                link_line = LINK_line( line )
                self.link_records.append( link_line )
            
            if line[0:4] == "ATOM":
                # store the protein lines
                pdb_line = PDB_line( line )
                
                # unknown amino acid - skip the PDB
                if pdb_line.res_name() == "UNK":
                    unknown.append( line )
                    self.unknown_res_pdb_names.append( pdb_name )
                    break
                
                # skip PDBs that have deuterium as an element
                if pdb_line.element() == 'D':
                    deuterium.append( line )
                    self.deuterium_pdb_names.append( pdb_name )
                    break
                    
                # skip hydrogen atoms
                if pdb_line.element() != 'H':
                    # collect the residue's unique name
                    pro_res_name = pdb_line.res_name()
                    pro_res_chain = pdb_line.res_chain()
                    pro_res_num = str( pdb_line.res_num() )
                    uniq_pro_name = pro_res_name + '_' + pro_res_chain + '_' + pro_res_num
                    
                    # check occupancy level
                    if pdb_line.occupancy() != 1.00:
                        # skip this PDB if they don't use letters for the alternate location identifiers
                        if pdb_line.alt_loc() == '':
                            print "## Skipping", self.name, "because it did not specify a residue's alternate location code ##"
                            return False
                    
                        # keep the 'A' alternate location residues
                        if pdb_line.alt_loc() == 'A':
                            # instantiate a dictionary key for this specific protein residue using its unique name 
                            # will be filled with the non-hydrogen ATOM lines later
                            if uniq_pro_name not in self.protein.keys():
                                self.protein[ uniq_pro_name ] = []
                            self.protein_lines.append( pdb_line )
                            
                    # otherwise this protein residue is at full occupancy so store the line
                    else:
                        # instantiate a dictionary key for this specific protein residue using its unique name 
                        # will be filled with the non-hydrogen ATOM lines later
                        if uniq_pro_name not in self.protein.keys():
                            self.protein[ uniq_pro_name ] = []
                        self.protein_lines.append( pdb_line )

                
            if line[0:6] == "HETATM":
                # get each ligand residue name to see what it is and store it in the appropriate list
                pdb_line = PDB_line( line )
                lig_res_name = pdb_line.res_name()
                
                # check what each residue is by its name (water, metal, ligand)
                if lig_res_name == "HOH":
                    water.append( line )
                elif lig_res_name == "DOD":
                    water.append( line )
                elif lig_res_name in metal_list:
                    metals.append( line )
                elif lig_res_name in AA_list:
                    AA_lig.append( line )
                elif lig_res_name in nucleic_acid_list:
                    nuc_acid_lig.append( line )
                # unknown amino acid - skip the PDB
                elif lig_res_name == "UNK":
                    unknown.append( line )
                    self.unknown_res_pdb_names.append( pdb_name )
                    break
                # unknown nucleic acid - skip the PDB
                elif lig_res_name == 'N':
                    unknown.append( line )
                    self.unknown_res_pdb_names.append( pdb_name )
                    break
                else:
                    # skip PDBs that have deuterium as an element
                    if pdb_line.element() == 'D':
                        deuterium.append( line )
                        self.deuterium_pdb_names.append( pdb_name )
                        break
                    
                    # otherwise, keep going
                    if pdb_line.element() != 'H':
                        # get the ligand residue's unique name
                        lig_res_chain = pdb_line.res_chain()
                        lig_res_num = str( pdb_line.res_num() )
                        uniq_lig_name = lig_res_name + '_' + lig_res_chain + '_' + lig_res_num
                        
                        # if this is a modified amino acid residue
                        if uniq_lig_name in modres_protein_res_names:
                            # change the HETATM to ATOM in the line
                            line = pdb_line.line
                            line = line.replace( "HETATM", "ATOM  " )
                            pdb_line = PDB_line( line )
                            
                            # now this is a unique protein name, not ligand
                            uniq_pro_name = uniq_lig_name
                                
                            # check occupancy level
                            if pdb_line.occupancy() != 1.00:
                                # skip this PDB if they don't use letters for the alternate location identifiers
                                if pdb_line.alt_loc() == '':
                                    print "## Skipping", self.name, "because it did not specify a residue's alternate location code ##"
                                    return False
                                
                                # keep the 'A' alternate location residues
                                if pdb_line.alt_loc() == 'A':
                                    # instantiate a dictionary key for this specific protein residue using its unique name 
                                    # will be filled with the non-hydrogen ATOM lines later
                                    if uniq_pro_name not in self.protein.keys():
                                        self.protein[ uniq_pro_name ] = []
                                    self.protein_lines.append( pdb_line )
                                    
                            # otherwise this modified residue is at full occupancy
                            else:
                                # instantiate a dictionary key for this specific modified protein residue 
                                # will be filled with the non-hydrogen HETATM lines later
                                if uniq_pro_name not in self.protein.keys():
                                    self.protein[ uniq_pro_name ] = []
                                    
                                # append the altered PDB line to the protein lines
                                self.protein_lines.append( pdb_line )
                                
                        
                        # otherwise, this is a ligand residue
                        else:
                            # check occupancy level
                            if pdb_line.occupancy() != 1.00:
                                # skip this PDB if they don't use letters for the alternate location identifiers
                                if pdb_line.alt_loc() == '':
                                    print "## Skipping", self.name, "because it did not specify a residue's alternate location code ##"
                                    return False
                                
                                # keep the 'A' alternate location residues
                                if pdb_line.alt_loc() == 'A':
                                    # instantiate a dictionary key for this specific ligand residue
                                    # will be filled with the non-hydrogen HETATM lines later
                                    if uniq_lig_name not in self.ligand.keys():
                                        self.ligand[ uniq_lig_name ] = []
                                        
                                    # store the HETATM lines
                                    self.hetatm_lines.append( pdb_line )
                                    
                            # otherwise residue is at full occupancy so store the line
                            else:
                                # instantiate a dictionary key for this specific ligand residue
                                # will be filled with the non-hydrogen HETATM lines later
                                if uniq_lig_name not in self.ligand.keys():
                                    self.ligand[ uniq_lig_name ] = []
                                    
                                # store the HETATM lines
                                self.hetatm_lines.append( pdb_line )

                            
        # if there were unknown residues, skip
        if len( unknown ) != 0:
            print "## Skipping", self.name, "because it has unknown residues ##"
            return False
         
        # if there were PDBs with more than one model, skip
        if len( models ) != 0:
            print "## Skipping", self.name, "because it has more than one model ##"
            return False
 
        # if there were PDBs with more than deuterium, skip
        if len( deuterium ) != 0:
            print "## Skipping", self.name, "because it contains deuterium ##"
            return False
        
        # if there is no ligand, skip
        if len( self.hetatm_lines ) == 0:
            print "## Skipping", self.name, "because it does not have a ligand of interest ##"
            return False

        # move all ATOM and HETATM lines to their appropriate place in the dictionaries based on their unique names
        for pro_pdb_line in self.protein_lines:
            # get the residue's unique name again
            pro_res_name = pro_pdb_line.res_name()
            pro_res_chain = pro_pdb_line.res_chain()
            pro_res_num = str( pro_pdb_line.res_num() )
            uniq_pro_name = pro_res_name + '_' + pro_res_chain + '_' + pro_res_num
            
            # otherwise just a normal residue
            self.protein[ uniq_pro_name ].append( pro_pdb_line )
        

        for lig_pdb_line in self.hetatm_lines:
            lig_res_name = lig_pdb_line.res_name()
            lig_res_chain = lig_pdb_line.res_chain()
            lig_res_num = str( lig_pdb_line.res_num() )
            uniq_lig_name = lig_res_name + '_' + lig_res_chain + '_' + lig_res_num
            
            # otherwise just a normal residue
            self.ligand[ uniq_lig_name ].append( lig_pdb_line )

        # if user wants to ignore glycosylated proteins (proteins with a HETATM attached to them)
        if ignore_glycosylated_proteins: 
            self.covalently_bound_lig_residues = []
            
            # see if there are residues covalently bound to the protein
            self.covalently_bound_lig_residues = self.determine_covalently_bound_ligands( pdb_name, self.protein, self.ligand, self.link_records )
            
            # add PDB to list if it had a glycan
            if len( self.covalently_bound_lig_residues ) != 0:
                # add name of PDB to glycosylated_proteins list (to be dumped later)
                self.glycosylated_proteins.append( self.name )
            
                # remove covalently bound ligands from the list of unique ligand residue names
                for remove_this_lig in self.covalently_bound_lig_residues:
                    # using an if statement because a covalently bound ligand residue could be an amino acid
                    if remove_this_lig in self.ligand.keys():
                        self.ligand.pop( remove_this_lig )
        
                    # if there is no ligand after removing glycans, skip
                    if len( self.ligand.keys() ) == 0:
                        print "## Skipping", self.name, "because it did not have a ligand of interest after removing glycans ##"
                        return False
        
        # make a clean PDB file that can be kept at the user's request to see if the program is doing what they think it is
        # only need to make the clean file if the user wants to keep it
        if keep_clean_pdbs:
            # create and open a XXXX.clean.pdb file name in the pdb directory
            cur_dir = os.getcwd() + '/'
            pdb_dir = cur_dir + 'pdbs/'
            clean_pdb_filename = pdb_dir + pdb_name + ".clean.pdb"
            
            # write the protein and ligand lines that are kept after this round of splitting
            with open( clean_pdb_filename, 'wb' ) as pdb_fh:
                for atom_key in self.protein.keys():
                    for atom_line in self.protein[ atom_key ]:
                        pdb_fh.write( atom_line.line )
                for hetatm_key in self.ligand.keys():
                    for hetatm_line in self.ligand[ hetatm_key ]:
                        pdb_fh.write( hetatm_line.line )
        
        return True

    

    def get_ligand_residues( self, heavy_atoms, cutoff):
        # dictionary -- key: lig residue number, value: hetatm lines
        self.ligand_dict = {}
        
        # list of the 3-letter ligand names, not allowing repeats
        # used later to count the number of each ligand residue
        self.uniq_lig_res_names = []
        
        # list of only the 3-letter ligand names, allowing repeats - used for data analysis
        self.lig_res_names = []
        
        # ligand data holders
        self.num_ligand_residues = 0
        self.num_ligand_atoms = 0
        self.num_ligand_nonpolar_atoms = 0
        self.num_ligand_polar_atoms = 0
        self.num_ligand_metal_atoms = 0
        self.num_ligand_unk_atom_type = 0
        self.lig_num_nonpolar_atoms = {}
        self.lig_num_polar_atoms = {}
        self.lig_num_metal_atoms = {}
        self.lig_num_unk_atoms = {}
        
        # this makes a new dictionary of ligands only if there are more than the specified number of heavy atoms in the lig residue
        for uniq_lig in self.ligand.keys():
            if len( self.ligand[ uniq_lig ] ) >= heavy_atoms:
                self.ligand_dict[ uniq_lig ] = self.ligand[ uniq_lig ]
        
        # get the 3-letter names of each remaining ligand residue
        for uniq_lig in self.ligand_dict.keys():
            lig_res_name = uniq_lig[ 0:3 ]
            self.lig_res_names.append( lig_res_name )
            
            # add unique ligand residue names to the set
            if lig_res_name not in self.uniq_lig_res_names:
                self.uniq_lig_res_names.append( lig_res_name )
        
        # get number of ligand residues
        self.num_ligand_residues = len( self.ligand_dict.keys() )
        
        # count the number of heavy ligand atoms by looping through the dictionary
        # hydrogens were skipped so just count the number of values for each key
        for uniq_lig in self.ligand_dict:
            self.num_ligand_atoms += len( self.ligand_dict[ uniq_lig ] )
            
            # prepare the counter for nonpolar, polar, metal, and unknown atom types for each unique ligand residue
            self.lig_num_nonpolar_atoms[ uniq_lig ] = 0
            self.lig_num_polar_atoms[ uniq_lig ] = 0
            self.lig_num_metal_atoms[ uniq_lig ] = 0
            self.lig_num_unk_atoms[ uniq_lig ] = 0
            
            # count the number of nonpolar and polar atoms
            for pdb_line in self.ligand_dict[ uniq_lig ]:                
                # if the element is nonpolar
                if pdb_line.element() in nonpolar_atoms:
                    # total nonpolar lig atoms
                    self.num_ligand_nonpolar_atoms += 1
                    
                    # number of nonpolar atoms for this ligand residue
                    self.lig_num_nonpolar_atoms[ uniq_lig ] += 1
                    
                # if the element is polar
                elif pdb_line.element() in polar_atoms:
                    # total polar lig atoms
                    self.num_ligand_polar_atoms += 1
                    
                    # number of polar atoms for this ligand residue
                    self.lig_num_polar_atoms[ uniq_lig ] += 1
                    
                # if the element is metal
                elif pdb_line.element() in metal_list:
                    # total metal lig atoms
                    self.num_ligand_metal_atoms += 1
                    
                    # number of metal atoms for this ligand residue
                    self.lig_num_metal_atoms[ uniq_lig ] += 1
                    
                    
                else:
                    # this atom type is unknown - inform user
                    print "      * I didn't know what type of atom", "'%s'" %pdb_line.element(), "is. Please add it to the list"
                    # total unknown lig atoms
                    self.num_ligand_unk_atom_type += 1
                    
                    # number of unknown atoms for this ligand residue
                    self.lig_num_unk_atoms[ uniq_lig ] += 1
                    
        # stop if there were no ligand residues with the given heavy atom cutoff
        if self.num_ligand_residues == 0:
            print "## Skipping", self.name, "because it had no ligand residues left given the", heavy_atoms, "heavy atom cutoff ##"
            return False
        
        # otherwise this ligand passed user requirements
        else:
            print "  ", self.name,
            print "has", self.num_ligand_residues, 
            print "ligand residues with more than", heavy_atoms, "heavy atoms",
            print "and", self.num_ligand_atoms, "non-hydrogen atoms"
            
            # replace self.hetatm_lines with only the hetatm lines of the remaining ligand residues with the appropriate number of heavy atoms
            # self.hetatm_lines will be used elsewhere - so it's easier to keep this accurate with only the hetatm lines that will be used
            self.hetatm_lines = []
            for uniq_lig in self.ligand_dict.keys():
                for pdb_line in self.ligand_dict[ uniq_lig ]:
                    self.hetatm_lines.append( pdb_line.line )
            
            return True
        
        
        
    def calc_distance(self, vec1, vec2):
        # takes two lists of xyz coordinates of two atoms
        from math import sqrt, pow
        
        x1 = vec1[0]
        y1 = vec1[1]
        z1 = vec1[2]
        
        x2 = vec2[0]
        y2 = vec2[1]
        z2 = vec2[2]
        
        dist = sqrt( pow( x2 - x1, 2 ) + pow( y2 - y1, 2 ) + pow( z2 - z1, 2 ) )
        
        return dist
    
    
    
    def get_activesite( self, cutoff ):
        # overall activesite dictionary
        # key: unique protein name (resname_reschain_resnum), value: list of ATOM lines per residue
        self.activesite_dict = {}
        
        # unique ligand name (resname_reschain_resnum)
        #value: list of 3-letter amino acid names
        self.activesite_lig_pro_res_dict = {}
        
        # unique ligand name (resname_reschain_resnum)
        #value: list of atom_line for each AA ( to get atom count later )
        self.activesite_lig_pro_atms_dict = {}
        
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
        
        for uniq_lig_name in self.ligand_dict.keys():
            # list to store the 3 letter names of all of the protein residues within the cutoff distance of each ligand residue (by unique name)
            AS_names_in_activesite = []
            AS_atms_in_activesite = []
            
            for hetatm_line in self.ligand_dict[ uniq_lig_name ]:
                # extract coordinates
                x_lig = hetatm_line.res_x_coord()
                y_lig = hetatm_line.res_y_coord()
                z_lig = hetatm_line.res_z_coord()
                lig_xyz = [ x_lig, y_lig, z_lig ]
            
                for atom_line in self.protein_lines: 
                    pro_res_name = atom_line.res_name()
                    pro_res_chain = atom_line.res_chain()
                    pro_res_num = str( atom_line.res_num() )
                    pro_uniq_res = str( pro_res_name + '_' + pro_res_chain + '_' + pro_res_num )
                    
                    # extract coordinates
                    x_pro = atom_line.res_x_coord()
                    y_pro = atom_line.res_y_coord()
                    z_pro = atom_line.res_z_coord()
                    pro_xyz = [ x_pro, y_pro, z_pro ]
                    
                    # check the distance
                    if self.calc_distance( lig_xyz, pro_xyz ) <= cutoff:
                        # append the line if the unique protein residue has already been counted
                        if pro_uniq_res in self.activesite_residues:
                            if atom_line not in self.activesite_dict[ pro_uniq_res ]:
                                self.activesite_dict[ pro_uniq_res ].append( atom_line )
                            
                        # store all of the unique names of the protein residues within the activesite
                        # also, store all of the atom_lines for each unique protein residue in the activesite
                        else:
                            self.activesite_residues.append( pro_uniq_res )
                            self.activesite_dict[ pro_uniq_res ] = []
                            self.activesite_dict[ pro_uniq_res ].append( atom_line )
                            
                            # store the 3 letter name of the amino acid within the activesite
                            three_letter_name = pro_uniq_res[0:3]
                            AS_names_in_activesite.append( three_letter_name )
                            
                        # store atom_line for each unique amino acid within the activesite to get an atom count later
                        if atom_line not in AS_atms_in_activesite:
                            AS_atms_in_activesite.append( atom_line )
                        
            # store the list of the 3 letter names for the amino acid within the self.activesite_lig_pro_dict according to which ligand its near
            self.activesite_lig_pro_res_dict[ uniq_lig_name ] = AS_names_in_activesite
            self.activesite_lig_pro_atms_dict[ uniq_lig_name ] = AS_atms_in_activesite
                
        # get number of activesite residues
        self.num_activesite_res = len( self.activesite_dict.keys() )
        
        # count the number of activesite atoms
        for uniq_lig_name in self.activesite_lig_pro_atms_dict.keys():
            self.num_activesite_atms += len( self.activesite_lig_pro_atms_dict[ uniq_lig_name ] )
            
            # prepare the counter for nonpolar, polar, and unknown atom types for each unique activesite residue
            self.activesite_num_nonpolar_atoms[ uniq_lig_name ] = 0
            self.activesite_num_polar_atoms[ uniq_lig_name ] = 0
            self.activesite_num_unk_atoms[ uniq_lig_name ] = 0
            
            # count the number of nonpolar and polar atoms
            for pdb_line in self.activesite_lig_pro_atms_dict[ uniq_lig_name ]:
                # if element is nonpolar
                if pdb_line.element() in nonpolar_atoms:
                    # total nonpolar activesite atoms
                    self.num_activesite_nonpolar_atoms += 1
                    
                    # number of nonpolar atoms for this activesite residue
                    self.activesite_num_nonpolar_atoms[ uniq_lig_name ] += 1
                    
                # if element is polar
                elif pdb_line.element() in polar_atoms:
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
        

        
    def get_activesite_AA_composition(self):
        ALA = []
        CYS = []
        ASP = []
        GLU = []
        PHE = []
        GLY = []
        HIS = []
        ILE = []
        LYS = []
        LEU = []
        MET = []
        ASN = []
        PRO = []
        GLN = []
        ARG = []
        SER = []
        THR = []
        VAL = []
        TRP = []
        TYR = []
        
        # loop over each unique protein residue and get the first three characters and up the appropriate count
        for pro_res in self.activesite_residues:
            res_name = pro_res[0:3]
            if res_name == "ALA":
                ALA.append( res_name )
            if res_name == "CYS":
                CYS.append( res_name )
            if res_name == "ASP":
                ASP.append( res_name )
            if res_name == "GLU":
                GLU.append( res_name )
            if res_name == "PHE":
                PHE.append( res_name )
            if res_name == "GLY":
                GLY.append( res_name )
            if res_name == "HIS":
                HIS.append( res_name )
            if res_name == "ILE":
                ILE.append( res_name )
            if res_name == "LYS":
                LYS.append( res_name )
            if res_name == "LEU":
                LEU.append( res_name )
            if res_name == "MET":
                MET.append( res_name )
            if res_name == "ASN":
                ASN.append( res_name )
            if res_name == "PRO":
                PRO.append( res_name )
            if res_name == "GLN":
                GLN.append( res_name )
            if res_name == "ARG":
                ARG.append( res_name )
            if res_name == "SER":
                SER.append( res_name )
            if res_name == "THR":
                THR.append( res_name )
            if res_name == "VAL":
                VAL.append( res_name )
            if res_name == "TRP":
                TRP.append( res_name )
            if res_name == "TYR":
                TYR.append( res_name )
                
        # collect the number of polar and nonpolar atoms
        percentage_activesite_nonpolar = round( float( self.num_activesite_nonpolar_atoms ) / float( self.num_activesite_atms ), 3 )
        percentage_activesite_polar = round( float( self.num_activesite_polar_atoms ) / float( self.num_activesite_atms ), 3 )
        percentage_ligand_nonpolar = round( float( self.num_ligand_nonpolar_atoms ) / float( self.num_ligand_atoms ), 3 )
        percentage_ligand_polar = round( float( self.num_ligand_polar_atoms ) / float( self.num_ligand_atoms ), 3 )
        percentage_ligand_metal = round( float( self.num_ligand_metal_atoms ) / float( self.num_ligand_atoms ), 3 )

        # append all of the final data to the self.lists
        # because if the analysis got this far, that means there actually is data to collect
        self.AS_pdb_names.append( self.name )
        self.AS_lig_res.append( self.num_ligand_residues )
        self.AS_lig_atms.append( self.num_ligand_atoms )
        self.AS_activesite_res.append( self.num_activesite_res )
        self.AS_activesite_atms.append( self.num_activesite_atms )
        self.ALA.append( len( ALA ) )
        self.CYS.append( len( CYS ) )
        self.ASP.append( len( ASP ) )
        self.GLU.append( len( GLU ) )
        self.PHE.append( len( PHE ) )
        self.GLY.append( len( GLY ) )
        self.HIS.append( len( HIS ) )
        self.ILE.append( len( ILE ) )
        self.LYS.append( len( LYS ) )
        self.LEU.append( len( LEU ) )
        self.MET.append( len( MET ) )
        self.ASN.append( len( ASN ) )
        self.PRO.append( len( PRO ) )
        self.GLN.append( len( GLN ) )
        self.ARG.append( len( ARG ) )
        self.SER.append( len( SER ) )
        self.THR.append( len( THR ) )
        self.VAL.append( len( VAL ) )
        self.TRP.append( len( TRP ) )
        self.TYR.append( len( TYR ) )
        self.AS_num_activesite_nonpolar_atoms.append( self.num_activesite_nonpolar_atoms )
        self.AS_num_activesite_polar_atoms.append( self.num_activesite_polar_atoms )
        self.AS_num_activesite_unk_atom_types.append( self.num_activesite_unk_atom_types )
        self.AS_num_ligand_nonpolar_atoms.append( self.num_ligand_nonpolar_atoms )
        self.AS_num_ligand_polar_atoms.append( self.num_ligand_polar_atoms )
        self.AS_num_ligand_metal_atoms.append( self.num_ligand_metal_atoms )
        self.AS_num_ligand_unk_atom_types.append( self.num_ligand_unk_atom_type )



    def get_activesite_AA_composition_per_lig_res(self):
        # goes through each unique ligand residue and counts the number of each amino acid within the cutoff distance around it
        for uniq_lig_name in self.activesite_lig_pro_res_dict.keys():
            # append the pdb names to the data list
            self.AS_pdb_names_per_lig.append( self.name )
            
            # append information about each ligand residue
            self.AS_lig_res_names_per_lig.append( uniq_lig_name.split( '_' )[0] )
            self.AS_lig_uniq_res_names_per_lig.append( uniq_lig_name )
            self.AS_lig_atms_per_lig.append( len( self.ligand_dict[ uniq_lig_name ] ) )
            
            # count and append the number of nonpolar and polar ligand atoms
            num_nonpolar_lig_atoms = self.lig_num_nonpolar_atoms[ uniq_lig_name ]
            num_polar_lig_atoms = self.lig_num_polar_atoms[ uniq_lig_name ]
            num_metal_lig_atoms = self.lig_num_metal_atoms[ uniq_lig_name ]
            num_unk_lig_atoms = self.lig_num_unk_atoms[ uniq_lig_name ]
                    
            self.AS_lig_num_ligand_nonpolar_atoms.append( num_nonpolar_lig_atoms )
            self.AS_lig_num_ligand_polar_atoms.append( num_polar_lig_atoms )
            self.AS_lig_num_ligand_metal_atoms.append( num_metal_lig_atoms )
            self.AS_lig_num_ligand_unk_atoms.append( num_unk_lig_atoms )
            
            # append information about all the activesite residues
            self.AS_activesite_res_per_lig.append( len( self.activesite_lig_pro_res_dict[ uniq_lig_name ] ) )
            self.AS_activesite_atms_per_lig.append( len( self.activesite_lig_pro_atms_dict[ uniq_lig_name ] ) )
            
            # count and append the number of nonpolar and polar activesite atoms
            num_nonpolar_activesite_atoms = self.activesite_num_nonpolar_atoms[ uniq_lig_name ]
            num_polar_activesite_atoms = self.activesite_num_polar_atoms[ uniq_lig_name ]
            num_unk_activesite_atoms = self.activesite_num_unk_atoms[ uniq_lig_name ]
            
            self.AS_lig_num_activesite_nonpolar_atoms.append( num_nonpolar_activesite_atoms )
            self.AS_lig_num_activesite_polar_atoms.append( num_polar_activesite_atoms )
            self.AS_lig_num_activesite_unk_atoms.append( num_unk_activesite_atoms )
            
            # count the number of amino acid residues around the ligand and append to data lists
            self.ALA_per_lig.append( self.activesite_lig_pro_res_dict[ uniq_lig_name ].count( "ALA" ) )
            self.CYS_per_lig.append( self.activesite_lig_pro_res_dict[ uniq_lig_name ].count( "CYS" ) )
            self.ASP_per_lig.append( self.activesite_lig_pro_res_dict[ uniq_lig_name ].count( "ASP" ) )
            self.GLU_per_lig.append( self.activesite_lig_pro_res_dict[ uniq_lig_name ].count( "GLU" ) )
            self.PHE_per_lig.append( self.activesite_lig_pro_res_dict[ uniq_lig_name ].count( "PHE" ) )
            self.GLY_per_lig.append( self.activesite_lig_pro_res_dict[ uniq_lig_name ].count( "GLY" ) )
            self.HIS_per_lig.append( self.activesite_lig_pro_res_dict[ uniq_lig_name ].count( "HIS" ) )
            self.ILE_per_lig.append( self.activesite_lig_pro_res_dict[ uniq_lig_name ].count( "ILE" ) )
            self.LYS_per_lig.append( self.activesite_lig_pro_res_dict[ uniq_lig_name ].count( "LYS" ) )
            self.LEU_per_lig.append( self.activesite_lig_pro_res_dict[ uniq_lig_name ].count( "LEU" ) )
            self.MET_per_lig.append( self.activesite_lig_pro_res_dict[ uniq_lig_name ].count( "MET" ) )
            self.ASN_per_lig.append( self.activesite_lig_pro_res_dict[ uniq_lig_name ].count( "ASN" ) )
            self.PRO_per_lig.append( self.activesite_lig_pro_res_dict[ uniq_lig_name ].count( "PRO" ) )
            self.GLN_per_lig.append( self.activesite_lig_pro_res_dict[ uniq_lig_name ].count( "GLN" ) )
            self.ARG_per_lig.append( self.activesite_lig_pro_res_dict[ uniq_lig_name ].count( "ARG" ) )
            self.SER_per_lig.append( self.activesite_lig_pro_res_dict[ uniq_lig_name ].count( "SER" ) )
            self.THR_per_lig.append( self.activesite_lig_pro_res_dict[ uniq_lig_name ].count( "THR" ) )
            self.VAL_per_lig.append( self.activesite_lig_pro_res_dict[ uniq_lig_name ].count( "VAL" ) )
            self.TRP_per_lig.append( self.activesite_lig_pro_res_dict[ uniq_lig_name ].count( "TRP" ) )
            self.TYR_per_lig.append( self.activesite_lig_pro_res_dict[ uniq_lig_name ].count( "TYR" ) )
            
        
        
        
    def count_contacts( self, cutoff ):
        # must have already found all residues within the activesite
        # ligand to protein!!
        self.polar_polar = 0
        self.polar_nonpolar = 0
        self.nonpolar_polar = 0
        self.nonpolar_nonpolar = 0
        self.unk_contact = 0
        # ligxyz_proxyz xyz coordinates unique for every atom, best way to collect unique contacts made
        self.uniq_contact_list = []
        
        for uniq_lig_name in self.ligand_dict.keys():
            # add pdb name and ligand name to list
            self.CC_per_lig_pdb_names.append( self.name )
            self.CC_per_lig_lig_names.append( uniq_lig_name )
            self.CC_per_lig_lig_atms.append( len( self.ligand_dict[ uniq_lig_name ] ) )
            self.CC_per_lig_activesite_atms.append( len( self.activesite_lig_pro_atms_dict[ uniq_lig_name ] ) )
            
            # make empty counters
            lig_polar_polar = 0
            lig_polar_nonpolar = 0
            lig_nonpolar_polar = 0
            lig_nonpolar_nonpolar = 0
            lig_unk_contact = 0
            
            # for every ligand residue
            for lig_pdb_line in self.ligand_dict[ uniq_lig_name ]:
                # get ligand atom xyz
                x_lig = lig_pdb_line.res_x_coord()
                y_lig = lig_pdb_line.res_y_coord()
                z_lig = lig_pdb_line.res_z_coord()
                lig_xyz = [ x_lig, y_lig, z_lig ]
                lig_xyz_str = str( x_lig ) + '_' + str( y_lig ) + '_' + str( z_lig )
                
                for pro_pdb_line in self.activesite_lig_pro_atms_dict[ uniq_lig_name ]:
                    # get protein atom xyz
                    x_pro = pro_pdb_line.res_x_coord()
                    y_pro = pro_pdb_line.res_y_coord()
                    z_pro = pro_pdb_line.res_z_coord()
                    pro_xyz = [ x_pro, y_pro, z_pro ]
                    pro_xyz_str = str( x_pro ) + '_' + str( y_pro ) + '_' + str( z_pro )
                    
                    # check atomic distance
                    if self.calc_distance( lig_xyz, pro_xyz ) < cutoff:
                        # check to see that this contact has not yet been counted
                       uniq_contact = lig_xyz_str + '.' + pro_xyz_str
                       if uniq_contact not in self.uniq_contact_list:
                           # append unique contact to list
                           self.uniq_contact_list.append( uniq_contact )
                           
                           # check contact type
                           if lig_pdb_line.element() in polar_atoms or lig_pdb_line.element() in metal_list:
                               # polar polar
                               if pro_pdb_line.element() in polar_atoms or pro_pdb_line.element() in metal_list:
                                   lig_polar_polar += 1
                                   
                               # polar nonpolar
                               elif pro_pdb_line.element() in nonpolar_atoms:
                                   lig_polar_nonpolar += 1
                                   
                               # unknown
                               else:
                                   lig_unk_contact += 1
                                   
                           elif lig_pdb_line.element() in nonpolar_atoms:
                               # nonpolar polar
                               if pro_pdb_line.element() in polar_atoms or pro_pdb_line.element() in metal_list:
                                   lig_nonpolar_polar += 1
                                   
                               # nonpolar nonpolar
                               elif pro_pdb_line.element() in nonpolar_atoms:
                                   lig_nonpolar_nonpolar += 1
                               
                               # unknown
                               else:
                                   lig_unk_contact += 1
                                   
                           # both unknown
                           else:
                               lig_unk_contact += 1
                            
            # append this info the per ligand data frame lists
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
        self.CC_pdb_names.append( self.name[0:4] )
        self.CC_lig_atms.append( self.num_ligand_atoms )
        self.CC_activesite_atms.append( self.num_activesite_atms )
        self.CC_pp_contacts.append( self.polar_polar )
        self.CC_pn_contacts.append( self.polar_nonpolar )
        self.CC_np_contacts.append( self.nonpolar_polar )
        self.CC_nn_contacts.append( self.nonpolar_nonpolar )
        self.CC_unk_contacts.append( self.unk_contact )
        



######################
#### RUNS PROGRAM ####
######################

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Use Python to count contacts.")
    parser.add_argument("pdb_name_list", help="a file of the pdbs to be analyzed")
    parser.add_argument("--ignore_glycosylated_proteins", "-i", action="store_true", help="do you want to skip PDBs that have a covalently attached HETATM group? This is most likely a glycan")
    parser.add_argument("--cutoff", "-c", type=int, default=5, help="how big do you want the activesite cutoff to be, in angstroms? default = 5")
    parser.add_argument("--heavy_atoms", "-ha", type=int, default=10, help="how many heavy atoms does a HETATM residue need to be considered a ligand? default = 10")
    parser.add_argument("--download_pdbs", "-d", action="store_true", help="do you need to download the pdbs from the database?")
    parser.add_argument("--keep_cifs", action="store_true", help="do you want to keep the cif files you download?")
    parser.add_argument("--keep_pdbs", action="store_true", help="do you want to keep the pdbs you download?")
    parser.add_argument("--keep_clean_pdbs", action="store_true", help="do you want to keep the cleaned-up version of the pdbs you are working with?")
    input_args = parser.parse_args()
