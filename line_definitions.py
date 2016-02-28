#!/usr/bin/python
__author__ = "morganlnance"


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
