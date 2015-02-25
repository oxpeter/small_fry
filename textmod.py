#!/usr/bin/env python

""" testing the ability to modify files"""

import re
import sys

from tempfile import mkstemp
from shutil import move
from os import remove, close


def window_wiping(file_path):
    #Create temp file
    fh, abs_path = mkstemp()
    new_file = open(abs_path,'w')
    old_file = open(file_path)
    
    for line in old_file:
        new_file.write(line.replace("\r", "\n"))
    
    #close temp file
    new_file.close()
    close(fh)
    old_file.close()
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)
    


def addname(file_path, accession_no, gene_name):
    #Create temp file
    fh, abs_path = mkstemp()
    new_file = open(abs_path,'w')
    old_file = open(file_path)
    
    for line in old_file:
        patt = "ID=" + accession_no
        geneid = re.search(patt, line)    
        
        if geneid is not None:
            oldline = line.strip()
            newline = oldline + "Name=" + gene_name + ";\n"
            new_file.write(newline)
        else:
            new_file.write(line)
    
    
    #close temp file
    new_file.close()
    close(fh)
    old_file.close()
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)
    
    
if __name__ == '__main__':
    cmmd, filename = sys.argv
    window_wiping(filename)