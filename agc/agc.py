#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""
#from operator import itemgetter
#import random
import gzip
import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Jules Collat"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Jules COllat"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Jules Collat"
__email__ = "collatjule@eisti.eu"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path

def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest = 'amplicon_file', type = isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest = 'minseqlen', type = int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest = 'mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest = 'chunk_size', type = int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest = 'kmer_size', type = int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest = 'output_file', type = str,
                        default="OTU.fasta", help = "Output file")
    return parser.parse_args()

def read_fasta(amplicon_file, minseqlen):
    """
     retourne un générateur de séquences de longueur I >= minseqlen
    """
    seq = ''
    with gzip.open(amplicon_file, "rt") as ampliconfile :
        for line in ampliconfile:
            line = line.replace("\n","")
            if   line != '' and line[:1] != ">" :
                seq+=line
            else:
                if len(seq) >= minseqlen:
                    yield seq
                seq =''   
                
def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    """
     retourne un générateur des séquences uniques ayant une occurrence O >= mincount ainsi que leur occurrence
    """
    dict_seq = {}
    generator = read_fasta(amplicon_file,minseqlen)
    for sequence in generator:
        if sequence not in dict_seq:
            dict_seq[sequence] = 1
        else:
            dict_seq[sequence]+= 1
    
    for seq, count in sorted(dict_seq.items(), key=lambda item: item[1], reverse=True):
        if count >= mincount:
            yield[seq,count]
    
def get_chunks(sequence, chunk_size):
    """
     retourne une liste de sous-séquences de taille I non chevauchantes.
    """
    liste_segment = []
    for i in range(0, len(sequence), chunk_size):
        if i+chunk_size <= len(sequence):
            liste_segment.append(sequence[i:i+chunk_size])
    if len(liste_segment)>=4:
        return liste_segment
    pass

def get_unique(ids):
    return {}.fromkeys(ids).keys()

def common(lst1, lst2): 
    return list(set(lst1) & set(lst2))

def cut_kmer(sequence, kmer_size):
    """
     retourne un générateur de tous les mots de longueur k présents dans cette séquence.
    """
    for i in range (len(sequence) - kmer_size +1):
        yield sequence[i :i+kmer_size]
    
def get_identity(alignment_list):
    """
    retourne le pourcentage d'identité entre les deux séquences.
    """
    nuc_identique = 0
    for i in range(len(alignment_list[0])):
        if alignment_list[0][i] == alignment_list[1][i]:
            nuc_identique += 1
    
    return (nuc_identique/len(alignment_list[0])) * 100

def get_unique_kmer(kmer_dict : dict, sequence : str, id_seq : int, kmer_size : int):   
    """
    retourne un dictionnaire de kmer contenant les kmers uniques présents dans chaque séquence pour 
    une longueur de donnée de kmer
    """
    for kmer in cut_kmer(sequence, kmer_size):
        if kmer not in kmer_dict:
            kmer_dict[kmer] = []
        kmer_dict[kmer].append(id_seq)
    return kmer_dict
    
def search_mates(kmer_dict : dict, sequence : str, kmer_size : int):
    """
    retourne les 8 séquences les plus similaires à notre séquence entrée
    """
    return [i[0] for i in Counter([ids for kmer in cut_kmer(sequence, kmer_size)
        if kmer in kmer_dict for ids in kmer_dict[kmer]]).most_common(8)]
    

def detect_chimera(perc_identity_matrix):
    """
    retourne un booléen indiquant si la séquence candidate est une chimère(True)
    ou ne l'est pas(False)
    """
    pass

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """
    retourne un générateur des séquences non chimérique au format : yield[sequence, count]
    """
    pass

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """
    retourne une liste d'OTU, cette liste indiquera pour chaque séquence son occurrence(count)
    """
    pass
    
def fill(text, width=80):
    """
    Split text with a line return to respect fasta format
    """
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
    """
     affiche les OTU au format 
     >OTU_{numéro partant de 1} occurrence:{nombre d'occurrence à la déréplication}
     {séquence au format fasta}
    """
    with open(output_file, "w") as outputfile:
        for i,OTU in enumerate(OTU_list):
            outputfile.write(">OTU_{} occurrence:{}\n".format(i+1, OTU[1]))
            outputfile.write("{}\n".format(fill(OTU[0])))
            
    
#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()


if __name__ == '__main__':
    main()