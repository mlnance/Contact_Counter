#!/usr/bin/python

from gzip import GzipFile
from StringIO import StringIO
import os
import urllib2
from io import BytesIO
from os.path import isfile
try:
    from colorama import Fore, Style
except:
    pass



def download_pdb(pdb_id, dest_dir):
    '''
    downloads pdbs, from yifans rosettacm py script
    :param pdb_id: 4letter accession code
    :param dest_dir:
    :return:
    '''
    try:
        pdb_name= pdb_id.split('_')[0]
    except:
        pass

    dest = '%s.pdb.gz' % ( pdb_name.lower() )
    filename = dest_dir + "/" + dest
    if not isfile(filename[:-3]):
        url = 'http://www.rcsb.org/pdb/files/%s.pdb.gz' % ( pdb_name.upper() )
        #curl_cmd = 'curl -o %s/%s %s' % ( dest_dir, dest, url )
        wget_cmd = 'wget --quiet %s -O %s' % ( url, filename )
        #try:
        #    log_lines = os.popen( curl_cmd ).readlines()
        #except:
        log_lines = os.popen( wget_cmd ).readlines()

        filename = dest_dir + "/" + dest
        # print filename
        os.system("gunzip -f %s" %filename)

        if ( os.path.isfile( filename[:-3] )):
            return dest[:-3]
        else:
            try:
                text = "~~Error: cannot download PDB %s!" %pdb_name
                print(Fore.BLUE + text + Style.RESET_ALL)
            except:
                print "~~Error: cannot download PDB %s!" %pdb_name
                

    return dest[:-3]



def download_cif_file(pdb_id, dest_dir):
    '''
    downloads cif file, from yifans rosettacm py script
    :param pdb_id: 4letter accession code
    :param dest_dir:
    :return:
    '''
    try:
        pdb_name= pdb_id.split('_')[0]
    except:
        pass

    dest = '%s.cif.gz' % ( pdb_name.lower() )
    filename = dest_dir + "/" + dest
    if not isfile(filename[:-3]):
        url = 'http://www.rcsb.org/pdb/files/%s.cif.gz' % ( pdb_name.upper() )
        #curl_cmd = 'curl -o %s/%s %s' % ( dest_dir, dest, url )
        wget_cmd = 'wget --quiet %s -O %s' % ( url, filename )
        #try:
        #    log_lines = os.popen( curl_cmd ).readlines()
        #except:
        log_lines = os.popen( wget_cmd ).readlines()

        filename = dest_dir + "/" + dest
        # print filename
        os.system("gunzip -f %s" %filename)

        if ( os.path.isfile( filename[:-3] )):
            return dest[:-3]
        else:
            print "Error: cannot download PDB %s!"%pdb_name

    return dest[:-3]
