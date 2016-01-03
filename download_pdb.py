#!/usr/bin/python

from gzip import GzipFile
from StringIO import StringIO
import os
import urllib2
from io import BytesIO
from os.path import isfile
import argparse

parser = argparse.ArgumentParser(description="Use Rosetta to make point mutations.")
parser.add_argument("pdb_id", help="the filename of the PDB structure to serve as the native structure")
parser.add_argument("dest_dir", help="where do you want your mutants to be dumped?")
input_args = parser.parse_args()

try:
    pdb_name = input_args.pdb_id.split('_')[0]
except:
    pass

dest = '%s.pdb.gz' % ( pdb_name.lower() )
filename = input_args.dest_dir + "/" + dest
if not isfile(filename[:-3]):
    url = 'http://www.rcsb.org/pdb/files/%s.pdb.gz' % ( pdb_name.upper() )
    #curl_cmd = 'curl -o %s/%s %s' % ( input_args.dest_dir, dest, url )
    wget_cmd = 'wget --quiet %s -O %s' % ( url, filename )
    #try:
    #    log_lines = os.popen( curl_cmd ).readlines()
    #except:
    log_lines = os.popen( wget_cmd ).readlines()
    
    filename = input_args.dest_dir + "/" + dest
    # print filename
    os.system("gunzip -f %s" %filename)
    
    if not ( os.path.isfile( filename[:-3] )):
        print "Error: cannot download PDB %s!"%pdb_name
