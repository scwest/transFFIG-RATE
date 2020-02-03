'''
Sean West 
31 January 2020

Simple code to do a single MSA run on a single gene.
'''

'''
Sean West 
31 January 2020

Simple code to do a single MSA run on a single gene.
'''

import sys 
import subprocess
from Bio import Phylo
import collections

class Transffig_mafft():
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
        command = ['mafft', '--distout', temp_fasta_filename, '>', temp_tree_filename]
        subprocess.call(command)
        
        # process mafft out into a distance matrix
        with open(temp_tree_filename, 'r') as infile,\
             open(distance_filename, 'w') as outfile:
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
        return

def smain():
    stick = Transffig_mafft()
    stick.main()
    return

