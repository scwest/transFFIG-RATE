'''
Sean West 
31 January 2020

Simple code to do a single MSA run on a single gene.

For information on MAFFT, see:
    Katoh, K., & Standley, D. M. (2013). MAFFT multiple 
    sequence alignment software version 7: improvements 
    in performance and usability. Molecular biology and 
    evolution, 30(4), 772-780.
'''

import sys 
import subprocess
from Bio import Phylo
import collections

class Transffig_mafft():
    def main(self):
        fasta_filename, distance_filename, unique_number = sys.argv[1:]
        
        # names of all involved temporary files (to be removed by end of script)
        temp_fasta_filename = 'mafft_fasta_temp{}.fa'.format(unique_number)
        temp_tree_filename = 'mafft_temp_tree{}.txt'.format(unique_number)
        
        # preprocess FASTA so Phylo doesn't mess up later 
        # Phylo can not handle a FASTA title with more than one quantity
        # (remove all but the transcript name)
        with open(fasta_filename, 'r') as infile, open(temp_fasta_filename, 'w') as outfile:
            for line in infile:
                if line[0] == '>':
                    line = line.strip().split(' ')[0] + '\n'
                outfile.write(line)
        
        # run mafft        
        command = ['mafft', '--distout', temp_fasta_filename]
        with open(temp_tree_filename, 'w') as outfile:
            subprocess.call(command, stdout=outfile)
        
        # process mafft out into a distance matrix
        dmat = collections.defaultdict(dict)
        with open(temp_fasta_filename+'.hat2', 'r') as infile,\
             open(distance_filename.replace('.csv', '.mafft.csv'), 'w') as outfile:
            infile.readline(); infile.readline(); infile.readline()
            names = []
            tri = []
            for line in infile:
                if line[:2] == '  ':
                    names.append(line.strip().split('=')[1].split(' ')[0])
                elif line[0] == ' ':
                    tri += line.strip().split(' ')
            for i in range(len(names)-1):
                k = -(i+1)
                m = len(names) + k
                j = 0
                while j < m:
                    value = tri.pop(0)
                    dmat[names[k]][names[j]] = value
                    dmat[names[j]][names[k]] = value
                    j += 1 
            names = sorted(names)
            outfile.write(','.join(names) + '\n')
            for i in names:
                outfile.write(i)
                for j in names:
                    if i == j:
                        outfile.write(',0')
                    else:
                        outfile.write(',' + dmat[i][j])
                outfile.write('\n')

        # remove any temporary files
        subprocess.call(['rm', temp_fasta_filename])
        subprocess.call(['rm', temp_tree_filename])
        subprocess.call(['rm', temp_fasta_filename+'.hat2'])
        return

def smain():
    '''
    This section is so that the transffig_mafft command can be run directly,
    without invoking transffig_rate or transffig_prep.
    
    This is also the reason that transffig_mafft is in its own package.
    This is useful when errors are produced for individual MSA jobs during
    transffig_prep.
    '''
    stick = Transffig_mafft()
    try:
        stick.main()
    except:
        print('Errors in transffig_mafft.')
    return

