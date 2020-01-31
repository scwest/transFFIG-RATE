'''
Sean West 
31 January 2020

Simple code to do a single MSA run on a single gene.
'''

import sys 
import subprocess
from Bio import Phylo
import collections

class Transffig_tcoffee():
    def main(self):
        fasta_filename, distance_filename, unique_number = sys.argv[1:]
        
        # names of all involved temporary files (to be removed by end of script)
        temp_fasta_filename = 'fasta_temp{}.fa'.format(unique_number)
        temp_tree_filename = 'temp_tree{}.txt'.format(unique_number)
        
        # preprocess fasta so Phylo doesn't mess up later 
        ## (remove all but the transcript name)
        with open(fasta_filename, 'r') as infile, open(temp_fasta_filename, 'w') as outfile:
            for line in infile:
                if line[0] == '>':
                    line = line.strip().split(' ')[0] + '\n'
                outfile.write(line)
        
        # run muscle
        command = ['t_coffee', temp_fasta_filename, '-newtree={}'.format(temp_tree_filename)]
        subprocess.call(command)
        
        # process t-coffee out into a distance matrix
        tree = Phylo.read(open(temp_tree_filename), 'newick')
        
        seqs = tree.get_terminals()
        dmat = collections.defaultdict(dict)
        for p1, p2 in combinations(seqs, 2):
            d = tree.distance(p1, p2)
            dmat[p1.name][p2.name] = d
            dmat[p2.name][p1.name] = d
            
        with open(distance_filename, 'w') as outfile:
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
        subprocess.call(['rm', temp_tree_filename])
        
        return
    
def smain():
    stick = Transffig_tcoffee()
    stick.main()



