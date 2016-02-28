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
sys.path.append( "utility" )
from download_files import download_pdb_file, download_cif_file
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

# list of residues to remove if they are a ligand
# ligand meaning it was defined as HETATM in the PDB
metal_list = [ 'B', 'K', 'V', 'Y', 'W', 'U', "LI", "BE", "NA", "MG", "AL", "SI", "CA", "SC", "TI", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "GA", "GE", "AS", "RB", "SR", "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN", "SB", "TE", "CS", "BA", "LU", "HF", "TA", "RE", "OS", "IR", "PT", "AU", "HG", "TI", "PB", "BI", "PO", "AT", "FR", "RA", "LR", "RF", "DB", "SG", "BH", "HS", "MT", "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", "AC", "TH", "PA", "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM", "MD", "NO" ]

# list of the three letter amino acid codes
AA_list = [ "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TRY", "VAL" ]

# list of the codes for the nucleic acids
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


class TER_line:
    def __init__(self, line):
        # only works for TER lines
        # line locations taken from the PDB databases's pdb format description page
        # ftp://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf
        # page 192
        self.line = line.rstrip( '/n' )
        
    def serial_number(self):
        return int( self.line[ 6:11 ].replace( ' ', '' ) )
        
    def res_name(self):
        return str( self.line[ 17:20 ].replace( ' ', '' ) )
        
    def res_chain(self):
        return str( self.line[ 21:22 ] )
        
    def res_seq(self):
        return int( self.line[ 22:26 ].replace( ' ', '' ) )
        
    def i_code(self):
        return str( self.line[ 26:27 ] )



class LINK_line:
    def __init__(self, line):
        # only works for LINK lines
        # line locations taken from the PDB databases's pdb format description page
        # ftp://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf
        # page 173
        self.line = line.rstrip( '/n' )
        
    def atom_name1(self):
        return str( self.line[ 12:16 ].replace( ' ', '' ) )
    
    def alt_loc1(self):
        return str( self.line[ 16:17].replace( ' ', '' ) )
    
    def res1_name(self):
        return str( self.line[17:20].replace( ' ', '' ) )
    
    def res1_chain(self):
        return str( self.line[21:22].replace( ' ', '' ) )
    
    def res1_seq(self):
        return int( self.line[22:26].replace( ' ', '' ) )
    
    def i_code1(self):
        return str( self.line[26:27].replace( ' ', '' ) )
    
    def atom_name2(self):
        return str( self.line[42:46].replace( ' ', '' ) )
    
    def alt_loc2(self):
        return str( self.line[46:47].replace( ' ', '' ) )
    
    def res2_name(self):
        return str( self.line[47:50].replace( ' ', '' ) )
    
    def res2_chain(self):
        return str( self.line[51:52].replace( ' ', '' ) )
    
    def res2_seq(self):
        return int( self.line[52:56].replace( ' ', '' ) )
    
    def i_code2(self):
        return str( self.line[56:57].replace( ' ', '' ) )
    
    def link_dist(self):
        return int( self.line[73:78].replace( ' ', '' ) )



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
            if connection_types[ii] == "covale" or connection_types[ii] == "metalc":
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



###########################
#### PDB CLEANING CODE ####
###########################

class Clean:
    def __init__( self ):
        self.working_dir = os.getcwd()
        
        # make on/off switch for skipping covalently bound ligand residues
        self.glycosylated_proteins = []
        
        # hold the name of PDBs that contain undesirables like multiple models or an unknown HETATM
        self.unknown_res_pdb_names = []
        self.multiple_models_pdb_names = []
        self.deuterium_pdb_names = []
        
        # hold unique ligand residue names for post-analysis
        self.uniq_lig_three_letter_codes_kept = []
        self.uniq_lig_three_letter_codes_skipped_size = []
        self.uniq_lig_three_letter_codes_skipped_link = []
        
        

    def instantiate_pdb_info_holders( self ):
        # instantiate lists that will hold relevant PDB data
        self.protein_lines = []
        self.hetatm_lines = []
        self.ter_lines = []
        self.link_records = []
        self.protein = {}
        self.ligand = {}        
        self.covalently_bound_lig_residues = []
        
        # dictionary -- key: unique ligand name, value: hetatm lines
        self.ligands_with_heavy_atom_cutoff = {}
        
        # ligand data holders for print PDB information to screen
        self.num_ligand_residues = 0
        self.num_ligand_atoms = 0
        
        return True



    def download_pdb( self, pdb_name ):
        """
        Uses a utility function to download a PDB from the internet based on the four letter accession code
        :param pdb_name: str( four-letter PDB code )
        :return: str( full path name to the PDB file )
        """
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
                pdb = download_pdb_file( pdb_name[:4], pdb_dir )
            except:
                pass
            
            # return the full path name to the PDB file
            return pdb_dir + pdb
    
        return True
    
    
    
    def download_cif( self, pdb_name ):
        """
        Uses a utility function to download the cif file of a PDB from the internet based on the four letter accession code
        :param pdb_name: str( four-letter PDB code )
        :return: str( full path name to the .cif file )
        """
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
        
        return True



    def pymol_clean( self, pdb_filename ):
        """
        Uses PyMOL to "clean" a PDB by removing waters and hydrogens.
        Opens up PyMOL via a command-line argument. Closes PyMOL upon completion
        :param pdb_filename:
        """
        # get the right names of the PDB files
        split_pdb_name = pdb_filename.split( '/' )[-1]
        pdb_name = split_pdb_name[:-4]
        
        # use pymol to remove water residues and remove hydrogens  -  ie. clean the file
        pymol_command = "%s -cqd 'load %s; remove resn hoh; remove resn dod; remove hydrogens %s; save %s'" %( pymol_exe_dir, pdb_filename, pdb_name, pdb_name )
        
        # run the PyMol command through the terminal
        os.popen( pymol_command )
        
        return True
    
    
    
    def get_uniq_connection_names_from_LINK_records( self, LINK_records ):
        unique_connection_names = []
        
        for link_line in LINK_records:
            res1_unique_name = link_line.res1_name() + '_' + link_line.res1_chain() + '_' + str( link_line.res1_seq() )
            res2_unique_name = link_line.res2_name() + '_' + link_line.res2_chain() + '_' + str( link_line.res2_seq() )
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
            
            # add each node according to whether its a HETATM or not
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
                        if g.node[ mynode ][ "HETATM" ]:
                            # if you haven't already added it
                            if mynode not in remove_these_ligs:
                                remove_these_ligs.append( mynode )
                                
        # return the unique names of the residues that are covalently linked to the protein
        return remove_these_ligs
    
    
    
    def determine_covalently_bound_ligands( self, pdb_filename, unique_protein_names, unique_hetatm_names, link_records ):
        # download the mmcif file
        cif_filename = self.download_cif( pdb_filename[:4] )
        
        # get _struct_conn lines to determine HETATM connections
        _struct_conn = cif_struct_conn_lines( cif_filename )
        
        # check to see if check_if_has_mmcif_dict() by using return value
        response = _struct_conn.check_if_has_mmcif_dict()
        
        # if this PDB has an mmcif file
        if response is True:
            # unique name = res1name_res1chain_res1num+res2name_res2chain_res2num
            unique_partner_names = _struct_conn.get_uniq_connection_names()
            print unique_partner_names
            
            # get list of ligand residues to remove
            remove_these_ligs = self.graph_out_residue_connections( unique_partner_names, unique_protein_names, unique_hetatm_names )
            
        # otherwise this PDB didn't have an mmcif file, use LINK records instead
        else:
            # unique name = res1name_res1chain_res1num+res2name_res2chain_res2num
            unique_partner_names = self.get_uniq_connection_names_from_LINK_records( link_records )
            
            # get list of ligand residues to remove
            remove_these_ligs = self.graph_out_residue_connections( unique_partner_names, unique_protein_names, unique_hetatm_names )
            
            
        return remove_these_ligs
    
    
    
    def split_pdb_file( self, pdb_filename, ignore_glycosylated_proteins ):
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
                
            if line[0:3] == "TER":
                # store TER lines to write to a clean PDB file
                ter_line = TER_line( line )
                self.ter_lines.append( ter_line )
            
            if line[0:4] == "ATOM":
                # store the protein lines
                pdb_line = PDB_line( line )
                
                # unknown amino acid - skip the PDB
                if pdb_line.res_name() == "UNK" or pdb_line.res_name() == "UNL":
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
                        # keep those without an identifier or with the code 'A'
                        if pdb_line.alt_loc() == '' or pdb_line.alt_loc() == 'A':
                            # instantiate a dictionary key for this specific protein residue using its unique name 
                            # will be filled with the non-hydrogen ATOM lines later
                            if uniq_pro_name not in self.protein.keys():
                                self.protein[ uniq_pro_name ] = []
                                
                            # store the protein lines
                            self.protein_lines.append( pdb_line )
                            
                    # otherwise this protein residue is at full occupancy so store the line
                    else:
                        # instantiate a dictionary key for this specific protein residue using its unique name 
                        # will be filled with the non-hydrogen ATOM lines later
                        if uniq_pro_name not in self.protein.keys():
                            self.protein[ uniq_pro_name ] = []
                            
                        # store the protein lines
                        self.protein_lines.append( pdb_line )
                        
                    
            if line[0:6] == "HETATM":
                # get each ligand residue name to see what it is and store it in the appropriate list
                pdb_line = PDB_line( line )
                lig_res_name = pdb_line.res_name()
                
                # check what each residue is by its name ( water, unknown, or a ligand etc )
                if lig_res_name == "HOH" or lig_res_name == "DOD":
                    water.append( line )
                elif lig_res_name in AA_list:
                    AA_lig.append( line )
                elif lig_res_name in nucleic_acid_list:
                    nuc_acid_lig.append( line )
                # unknown ligand - skip the PDB
                elif lig_res_name == "UNK" or lig_res_name == "UNL":
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
                                # keep those without an identifier or with the code 'A'
                                if pdb_line.alt_loc() == '' or pdb_line.alt_loc() == 'A':
                                    # instantiate a dictionary key for this specific protein residue using its unique name 
                                    # will be filled with the non-hydrogen ATOM lines later
                                    if uniq_lig_name not in self.ligand.keys():
                                        self.ligand[ uniq_lig_name ] = []
                                        
                                    # append the altered PDB line to the protein lines
                                    self.hetatm_lines.append( pdb_line )
                                    
                            # otherwise this modified residue is at full occupancy
                            else:
                                # instantiate a dictionary key for this specific modified protein residue 
                                # will be filled with the non-hydrogen HETATM lines later
                                if uniq_lig_name not in self.ligand.keys():
                                    self.ligand[ uniq_lig_name ] = []
                                    
                                # append the altered PDB line to the protein lines
                                self.hetatm_lines.append( pdb_line )
                            
                            
                        # otherwise, this is a ligand residue
                        else:
                            # check occupancy level
                            if pdb_line.occupancy() != 1.00:
                                # keep those without an identifier or with the code 'A'
                                if pdb_line.alt_loc() == '' or pdb_line.alt_loc() == 'A':
                                    # instantiate a dictionary key for this specific protein residue using its unique name 
                                    # will be filled with the non-hydrogen ATOM lines later
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
            print "## Skipping", pdb_name, "because it has unknown residues ##"
            return False
        
        # if there were PDBs with more than one model, skip
        if len( models ) != 0:
            print "## Skipping", pdb_name, "because it has more than one model ##"
            return False
        
        # if there were PDBs with more than deuterium, skip
        if len( deuterium ) != 0:
            print "## Skipping", pdb_name, "because it contains deuterium ##"
            return False
        
        # if there is no ligand, skip
        if len( self.hetatm_lines ) == 0:
            print "## Skipping", pdb_name, "because it does not have a ligand of interest ##"
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
            # see if there are residues covalently bound to the protein
            self.covalently_bound_lig_residues = self.determine_covalently_bound_ligands( pdb_name, self.protein, self.ligand, self.link_records )
            
            # add PDB to list if it had a glycan
            if len( self.covalently_bound_lig_residues ) != 0:
                # add name of PDB to glycosylated_proteins list (to be dumped later)
                self.glycosylated_proteins.append( pdb_name )
                
                # remove covalently bound ligands from the list of unique ligand residue names
                for remove_this_lig in self.covalently_bound_lig_residues:
                    # using an if statement because a covalently bound ligand residue could be an amino acid
                    if remove_this_lig in self.ligand.keys():
                        self.ligand.pop( remove_this_lig )
                        
                        # store the unique three letter code of the ligand residues that were skipped
                        three_letter_code = remove_this_lig.split( '_' )[0]
                        if three_letter_code not in self.uniq_lig_three_letter_codes_skipped_link:
                            self.uniq_lig_three_letter_codes_skipped_link.append( three_letter_code )
                        
                # if there is no ligand after removing glycans, skip
                if len( self.ligand.keys() ) == 0:
                    print "## Skipping", pdb_name, "because it did not have a ligand of interest after removing glycans ##"
                    return False
                
        return True
    
    
    
    def get_ligand_residues( self, heavy_atoms, pdb_filename ):
        # get the PDB name from the end of the full path given in pdb_filename
        split_pdb_name = pdb_filename.split( '/' )[-1]
        pdb_name = split_pdb_name[:-4].lower()
        
        # this makes a new dictionary of ligands only if there are more than the specified number of heavy atoms in the lig residue
        for uniq_lig in self.ligand.keys():
            # keep PDB based on its unique name and corresponding HETATM lines
            if len( self.ligand[ uniq_lig ] ) >= heavy_atoms:
                self.ligands_with_heavy_atom_cutoff[ uniq_lig ] = self.ligand[ uniq_lig ]
                
                # keep unique three letter codes for ligand residues kept
                three_letter_code = uniq_lig.split( '_' )[0]
                if three_letter_code not in self.uniq_lig_three_letter_codes_kept:
                    self.uniq_lig_three_letter_codes_kept.append( three_letter_code )
            else:                
                # otherwise it didn't pass the size cutoff
                # so keep unique three letter codes for ligand residues skipped 
                three_letter_code = uniq_lig.split( '_' )[0]
                if three_letter_code not in self.uniq_lig_three_letter_codes_skipped_size:
                    self.uniq_lig_three_letter_codes_skipped_size.append( three_letter_code )
               
        # get number of ligand residues
        self.num_ligand_residues = len( self.ligands_with_heavy_atom_cutoff.keys() )
        
        # count the number of heavy ligand atoms by looping through the dictionary
        # hydrogens were skipped so just count the number of values for each key
        for uniq_lig in self.ligands_with_heavy_atom_cutoff:
            self.num_ligand_atoms += len( self.ligands_with_heavy_atom_cutoff[ uniq_lig ] )
            
        # stop if there were no ligand residues with the given heavy atom cutoff
        if self.num_ligand_residues == 0:
            print "## Skipping", pdb_name, "because it had no ligand residues left given the", heavy_atoms, "heavy atom cutoff ##"
            return False
        
        # otherwise this ligand passed user requirements
        else:
            print "  ", pdb_name,
            print "has", self.num_ligand_residues, 
            print "ligand residues that each have equal to or more than", heavy_atoms, "heavy atoms",
            print "resulting in", self.num_ligand_atoms, "total heavy ligand atoms"
            
            # replace hetatm_lines with only the hetatm lines of the remaining ligand residues with the appropriate number of heavy atoms
            # hetatm_lines will be used elsewhere - so it's easier to keep this accurate with only the hetatm lines that will be used
            self.hetatm_lines = []
            for uniq_lig in self.ligands_with_heavy_atom_cutoff.keys():
                for pdb_line in self.ligands_with_heavy_atom_cutoff[ uniq_lig ]:
                    self.hetatm_lines.append( pdb_line.line )
                    
            return True
        
        return True
    
    
    
    def write_clean_pdb_file( self, pdb_filename ):
        # make a clean PDB file that can be kept at the user's request to see if the program is doing what they think it is
        # only need to make the clean file if the user wants to keep it
        # create and open a XXXX.clean.pdb file name in the pdb directory
        cur_dir = os.getcwd() + '/'
        pdb_dir = cur_dir + 'pdbs/'
        
        # get the four letter code in case a file path was given
        pdb_name = pdb_filename.split( '/' )[-1][:4]
        clean_pdb_filename = pdb_dir + pdb_name + ".clean.pdb"
        
        # put the ATOM, HETATM, and TER lines into a list to be sorted on atom number
        atom_num_pdb_lines = {}
        for atom_res_key in self.protein.keys():
            for atom_line in self.protein[ atom_res_key ]:
                atom_num_pdb_lines[ atom_line.atom_num() ] = atom_line.line
        for hetatm_res_key in self.ligand.keys():
            for hetatm_line in self.ligand[ hetatm_res_key ]:
                atom_num_pdb_lines[ hetatm_line.atom_num() ] = hetatm_line.line
        for ter_line in self.ter_lines:
            # serial_number == atom_number
            atom_num_pdb_lines[ ter_line.serial_number() ] = ter_line.line
                        
        # get the atom numbers in order from the dictionary
        ordered_atom_numbers = atom_num_pdb_lines.keys()
        ordered_atom_numbers.sort()
                
        # order the PDB lines based on their atom number
        ordered_pdb_lines = []
        for atom_num_key in ordered_atom_numbers:
            ordered_pdb_lines.append( atom_num_pdb_lines[ atom_num_key ] )
                    
        # write the protein and ligand lines that are kept after this round of splitting
        with open( clean_pdb_filename, 'wb' ) as pdb_fh:
            pdb_fh.writelines( ordered_pdb_lines )
                    
        return True
    
    
    
    
######################
#### RUNS PROGRAM ####
######################

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Use Python to count contacts.")
    parser.add_argument("pdb_name_list", help="a file of the pdbs to be analyzed")
    parser.add_argument("--ignore_glycosylated_proteins", "-i", action="store_true", help="do you want to skip PDBs that have a covalently attached HETATM group? This is most likely a glycan")
    parser.add_argument("--heavy_atoms", "-ha", type=int, default=10, help="how many heavy atoms does a HETATM residue need to be considered a ligand? default = 10")
    parser.add_argument("--download_pdbs", "-d", action="store_true", help="do you need to download the pdbs from the database?")
    parser.add_argument("--keep_cifs", action="store_true", help="do you want to keep the cif files you download?")
    parser.add_argument("--keep_pdbs", action="store_true", help="do you want to keep the pdbs you download?")
    parser.add_argument("--keep_clean_pdbs", action="store_true", help="do you want to keep the cleaned-up version of the pdbs you are working with?")
    input_args = parser.parse_args()
