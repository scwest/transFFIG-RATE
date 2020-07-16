'''
Sean West 
31 January 2020

Simple code to do a single MSA run on a single gene.

For information on Clustal Omega, see:
    Sievers, F., & Higgins, D. G. (2014). Clustal Omega, 
    accurate alignment of very large numbers of sequences. 
    In Multiple sequence alignment methods (pp. 105-116). 
    Humana Press, Totowa, NJ.

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
        
        # preprocess FASTA so Phylo doesn't mess up later 
        # Phylo can not handle a FASTA title with more than one quantity
        # (remove all but the transcript name)
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
             open(distance_filename.replace('.csv', '.clustalo.csv'), 'w') as outfile:
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
    '''
    This section is so that the clustalo command can be run directly,
    without invoking transffig_rate or transffig_prep.
    
    This is also the reason that transffig_clustalo is in its own package.
    This is useful when errors are produced for individual MSA jobs during
    transffig_prep.
    '''
    stick = Transffig_clustalo()
    try:
        stick.main()
    except:
        print('Errors in transffig_clustalo.')
    return


