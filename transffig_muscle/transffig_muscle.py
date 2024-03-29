'''
Sean West 
31 January 2020

Simple code to do a single MSA run on a single gene.

For information on MUSCLE, see:
    Edgar, R. C. (2004). MUSCLE: multiple sequence alignment 
    with high accuracy and high throughput. Nucleic acids 
    research, 32(5), 1792-1797.

'''

import sys 
import subprocess
from Bio import Phylo
import collections
from itertools import combinations

class Transffig_muscle():
    def main(self):
        fasta_filename, distance_filename, unique_number = sys.argv[1:]
        
        # names of all involved temporary files (to be removed by end of script)
        temp_fasta_filename = 'muscle_fasta_temp{}.fa'.format(unique_number)
        temp_empty_filename = 'muscle_temp_empty{}.txt'.format(unique_number)
        temp_tree_filename = 'muscle_temp_tree{}.txt'.format(unique_number)
        
        # preprocess FASTA so Phylo doesn't mess up later 
        # Phylo can not handle a FASTA title with more than one quantity
        # (remove all but the transcript n
        with open(fasta_filename, 'r') as infile, open(temp_fasta_filename, 'w') as outfile:
            for line in infile:
                if line[0] == '>':
                    line = line.strip().split(' ')[0] + '\n'
                outfile.write(line)
        
        # run muscle
        command = ['muscle', '-in', temp_fasta_filename, '-out', temp_empty_filename, \
                   '-tree2', temp_tree_filename]
        subprocess.call(command)
        
        # process muscle out into a distance matrix
        tree = Phylo.read(open(temp_tree_filename), 'newick')
        
        
        seqs = tree.get_terminals()
        dmat = collections.defaultdict(dict)
        for p1, p2 in combinations(seqs, 2):
            d = tree.distance(p1, p2)
            dmat[p1.name][p2.name] = d
            dmat[p2.name][p1.name] = d
            
        with open(distance_filename.replace('.csv', '.muscle.csv'), 'w') as outfile:
            names = sorted(dmat.keys())
            outfile.write(','.join(dmat) + '\n')
            for i in names:
                outfile.write(i)
                for j in names:
                    if i == j:
                        outfile.write(',0')
                    else:
                        outfile.write(',' + str(dmat[i][j]))
                outfile.write('\n')
        
        # remove any temporary files
        subprocess.call(['rm', temp_fasta_filename])
        subprocess.call(['rm', temp_fasta_filename.replace('.fa', '.aln')])
        subprocess.call(['rm', temp_fasta_filename.replace('.fa', 'html')])
        subprocess.call(['rm', temp_tree_filename])
        subprocess.call(['rm', temp_empty_filename])

        return
    
def smain():
    '''
    This section is so that the transffig_muscle command can be run directly,
    without invoking transffig_rate or transffig_prep.
    
    This is also the reason that transffig_muscle is in its own package.
    This is useful when errors are produced for individual MSA jobs during
    transffig_prep.
    '''
    stick = Transffig_muscle()
    try:
        stick.main()
    except:
        print('Errors in transffig_muscle.')
    return

