'''
Sean West 
31 January 2020

Simple code to do a single MSA run on a single gene.
'''

import sys 
import subprocess
from Bio import Phylo
import collections

class Transffig_clustalo():
    def main(self):
        fasta_filename, distance_filename, unique_number = sys.argv[1:]
        
        # names of all involved temporary files (to be removed by end of script)
        temp_fasta_filename = 'clustalo_fasta_temp{}.fa'.format(unique_number)
        temp_tree_filename = 'clustalo_temp_tree{}.txt'.format(unique_number)
        
        # preprocess fasta so Phylo doesn't mess up later 
        ## (remove all but the transcript name)
        with open(fasta_filename, 'r') as infile, open(temp_fasta_filename, 'w') as outfile:
            for line in infile:
                if line[0] == '>':
                    line = line.strip().split(' ')[0] + '\n'
                outfile.write(line)
        
        # run muscle
        command = ['clustalo', '--in={}'.format(temp_fasta_filename), \
                   '--distmat-out={}'.format(temp_tree_filename), '--force', '--full']
        subprocess.call(command)
        
        # process clustal omega out into a distance matrix
        with open(temp_tree_filename, 'r') as infile, \
             open(distance_filename.replace('.fa', '.clustalo.fa'), 'w') as outfile:
            infile.readline()
            names = []
            lines = []
            for line in infile:
                line = line.strip().split(' ', 1)
                names.append(line[0])
                lines.append((','.join(line[1].split(' '))).strip(',') + '\n')
            outfile.write(','.join(names) + '\n')
            for line in lines:
                n = names.pop(0)
                outfile.write(n + ',' + line)
        
        # remove any temporary files
        subprocess.call(['rm', temp_fasta_filename])
        subprocess.call(['rm', temp_tree_filename])
        return

def smain():
    stick = Transffig_clustalo()
    stick.main()
    return


