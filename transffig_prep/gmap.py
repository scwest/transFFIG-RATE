'''
Sean West
30 January 2020

This portion of the code handles:
 - processing GMAP output
 - preparing gene-specific fasta files
 - preparing a list of commands to run the MSAs

'''

import os
import collections
import sys
import subprocess
from copy import deepcopy

from transffig_prep import Gene
from transffig_prep import Gmap_chunk

class Gmap():
    def __init__(self):
        self.genes = {}
        self.current_gene_number = 1
        self.tran2sequence = {}
                        
    def absorb_hit(self, hit):
        '''
        New genes receive an 'arbitrary_<>' gene id.
        This is to emphasize that they are not real ids. That they may be reassigned to real
          gene ids if desired.
          
        ABSORB CHUNK DOES NOT ASSIGN SEQUENCES TO TRANSCRIPTS
        '''
        
        # Check to see all the genes that overlap the genomic location of the transcript
        overlapped_gene_names = []
        for gene_name, gene in self.genes.items():
            if self.isoverlap(gene, hit):
                overlapped_gene_names.append(gene_name)
            
        # Incorporate the transcript into existing (or new) genes
        ### If there is no gene that overlaps, we'll create a new gene
        ### create a new gene at this location
        if len(overlapped_gene_names) == 0:
            pass
            '''IGNORING TRANSCRIPTS NOT ASSIGNED TO A KNOWN GENE LOCATION
            # produce a new gene
            gene_name = 'arbitrary_{}'.format(self.current_gene_number)
            self.current_gene_number += 1
            t = collections.defaultdict(str)
            t[hit.name] = self.tran2sequence[hit.name]
            gene = Gene(name = gene_name,\
                        chromosome = hit.chromosome,\
                        strand = hit.strand,\
                        start = hit.start,\
                        end = hit.end,\
                        trans = t)
            self.genes[gene.name] = gene
            '''
        
        ### If we found 1 gene that the transcript overlaps;
        ### add to an existing gene
        elif len(overlapped_gene_names) == 1:
            gene_name = overlapped_gene_names.pop()
            self.genes[gene_name].start = min([self.genes[gene_name].start, hit.start])
            self.genes[gene_name].end = max([self.genes[gene_name].end, hit.end])
            self.genes[gene_name].trans[hit.name] = self.tran2sequence[hit.name]
        
        ### If we found multiple genes that the transcript overlaps;
        ### These are either genes that are on the same genomic location
        ### or they are the same gene. Either way, add them together.
        else:
            # for now, transcripts that overlap with reference genes will be represented
            #          and those that do not overlap with reference genes at all
            keeper_gene_names = [x for x in overlapped_gene_names if 'arbitrary' not in x]
            for gene_name in keeper_gene_names:
                self.genes[gene_name].start = min([self.genes[gene_name].start, hit.start])
                self.genes[gene_name].end = max([self.genes[gene_name].end, hit.end])
                self.genes[gene_name].trans[hit.name] = self.tran2sequence[hit.name]
            
            ''' IGNORING TRANSCRIPTS NOT ASSIGNED TO A KNOWN GENE LOCATION
            overlapped_gene_names = [x for x in overlapped_gene_names if 'arbitrary' in x]
            
            if len(overlapped_gene_names) == 1:
                gene_name = overlapped_gene_names.pop()
                self.genes[gene_name].start = min([self.genes[gene_name].start, hit.start])
                self.genes[gene_name].end = max([self.genes[gene_name].end, hit.end])
                self.genes[gene_name].trans[hit.name] = self.tran2sequence[hit.name]
            elif len(overlapped_gene_names) > 1:
                gene = Gene()
                gene.start = min([self.genes[x].start for x in overlapped_gene_names]+[hit.start])
                gene.end = max([self.genes[x].end for x in overlapped_gene_names]+[hit.end])
                gene.chromosome = hit.chromosome # they're all the same; only [0],[1] are guaranteed
                gene.strand = hit.strand
                
                #gene_names = [self.genes[x].name for x in overlapped_gene_names]
                keeper_gene_names = [x for x in overlapped_gene_names if 'arbitrary' not in x]
                if keeper_gene_names:
                    gene.name = '-'.join(keeper_gene_names)
                else:
                    gene.name = 'arbitrary_{}'.format(self.current_gene_number)
                    self.current_gene_number += 1
                
                for g in [self.genes[x] for x in overlapped_gene_names]:
                    for t, seq in g.trans.items():
                        gene.trans[t] = seq
                gene.trans[hit.name] = self.tran2sequence[hit.name]
                
                self.genes[gene.name] = gene
                
                for gene_name in overlapped_gene_names:
                    del self.genes[gene_name]
            '''
        return
        
    def check(self, storage_prefix):
        commands = []
        storage_filename = '{}commands.txt'.format(storage_prefix)
        
        if os.path.isfile(storage_filename):
            with open(storage_filename, 'r') as infile:
                for line in infile:
                    commands.append(line.strip().split(' '))
        return commands
    
    def command(self, gene, storage_prefix):
        commands = []
        for call in ['transffig_muscle', 'transffig_mafft', 'transffig_tcoffee', 'transffig_clustalo']:
            commands.append([call, gene.fa_filename, '{}distance_matrices/{}.csv'.format(storage_prefix, gene.name)])
        
        return commands
    
    def get_fasta_chunks(self, input_filename):
        with open(input_filename, 'r') as infile:
            chunk_of_text = infile.readline()
            for line in infile:
                if line[0] == '>':
                    yield chunk_of_text
                    chunk_of_text = line
                else:
                    chunk_of_text += line
        return chunk_of_text
    
    def get_gmap_chunks(self, input_filename):
        with open(input_filename, 'r') as infile:
            chunk_of_text = infile.readline()
            for line in infile:
                if line[0] == '>':
                    yield chunk_of_text
                    chunk_of_text = line
                else:
                    chunk_of_text += line
        return chunk_of_text
    
    def isoverlap(self, gene1, gene2):
        if gene1.chromosome == gene2.chromosome and gene1.strand == gene2.strand:
            if gene2.start > gene1.start and gene2.start < gene1.end:
                return True
            if gene2.end > gene1.start and gene2.end < gene1.end:
                return True
        return False
    
    def link_genes(self, reference_filename):
        '''
        The file must be setup as:
        <gene name>,<chromosome>,<start>,<end>,<strand>
        '''
        with open(reference_filename, 'r') as infile:
            c = 1
            for line in infile:
                sys.stdout.write('\r\t\tAdding Gene:\t{}'.format(c))
                sys.stdout.flush()
                c += 1
                
                line = line.strip().split(',')
                new_gene = Gene(name = line[0],\
                                chromosome = line[1],\
                                start = int(line[2]),\
                                end = int(line[3]),\
                                strand = int(line[4]))
                self.genes[new_gene.name] = deepcopy(new_gene)
                self.genes[new_gene.name].trans = collections.defaultdict(str)
            sys.stdout.write('\n')
        return
    
    def parse(self, storage_prefix, fasta_input_filename, gmap_output_filename, reference_filename):
        commands = []
        fasta_output_filename = '{}full.fa'.format(storage_prefix)
        
        # setup output directories
        fasta_dir = '{}gene_fastas'.format(storage_prefix)
        subprocess.call(['mkdir', fasta_dir])
        dm_dir = '{}distance_matrices'.format(storage_prefix)
        subprocess.call(['mkdir', dm_dir])
        
        # input all transcript sequences
        self.upload_transcript_fasta(fasta_input_filename)
        
        # check for an existing .fa that has the gmap already processed
        if os.path.isfile(fasta_output_filename):
            print('\tParsing Existing Commands.')
            self.parse_existing_fa(fasta_output_filename, storage_prefix)
        else:
            print('\tAdding initial set of reference genes.')
            if reference_filename:
                self.link_genes(reference_filename)
                
            print('\tLength of self.genes:\t{}'.format(len(self.genes)))
            
            
            print('\tParsing GMAP hits (chunks)')
            self.parse_gmap(gmap_output_filename, fasta_input_filename, fasta_output_filename)
            
            print('\tLength of self.genes:\t{}'.format(len(self.genes)))
                
            # write the genes to their individual fa's
            print('\tWriting gene-specific fasta files.')
            for gene_name, gene in self.genes.items():
                if len(gene.trans) < 3:
                    continue
                self.genes[gene_name].fa_filename = self.write_gene_fa(storage_prefix, gene)
                
            print('\tWriting new fasta with Genes to:\n\t\t{}'.format(fasta_output_filename))
            self.write_new_fasta(fasta_output_filename, fasta_input_filename)
            
        # write the commands for each gene
        i = 1
        t = len(self.genes)
        for gene in self.genes.values():
            if len(gene.trans) < 3:
                sys.stdout.write('\r\t\tGenes with < 3 transcripts: {} of {}'.format(i, t))
                sys.stdout.flush()
                i += 1
                continue
            #print(gene.name)
            #print(gene.trans)
            commands += self.command(gene, storage_prefix)
        print('')
        
        return commands 
    
    def add_gene_and_transcript(self, gene_name, tran_name, seq, storage_prefix):
        if gene_name in self.genes:
            self.genes[gene_name].trans[tran_name] = seq
        else:
            gene = self.Gene(name = gene_name)
            t = {tran_name:seq}
            gene.trans = t
            gene.fa_filename = '{}gene_fastas/{}.fa'.format(storage_prefix, gene_name)
            self.genes[gene_name] = gene
        return
    
    def parse_fa_chunk(self, chunk):
        lines = chunk.split('\n', 1)
        tran_name = lines[0].strip().split(' ')[0].replace('>', '')
        gene_name = lines[0].strip().split('gene:')[-1].split(' ')[0]
        print(lines)
        seq = lines[1]
        return tran_name, gene_name, seq
    
    def parse_existing_fa(self, fasta_output_filename, storage_prefix):
        with open(fasta_output_filename, 'r') as infile:
            chunk = infile.readline()
            for line in infile:
                if line[0] == '>':
                    tran_name, gene_name, seq = self.parse_fa_chunk(chunk)
                    self.add_gene_and_transcript(gene_name, tran_name, seq, storage_prefix)
                    chunk = line
                else:
                    chunk += line
            tran_name, gene_name, seq = self.parse_fa_chunk(chunk)
            self.add_gene_and_transcript(gene_name, tran_name, seq, storage_prefix)
        return
    
    def parse_gmap(self, gmap_output_filename, fasta_input_filename, fasta_output_filename):
        # get all the genes
        c = 1
        for chunk_of_text in self.get_gmap_chunks(gmap_output_filename):
            sys.stdout.write('\r\t\tChunk: {}'.format(c))
            sys.stdout.flush()
            c += 1
            
            chunk = Gmap_chunk(chunk_of_text)
            chunk.process()
            for hit in chunk.hits:
                self.absorb_hit(hit)
        print('')
                    
        return
    
    def upload_transcript_fasta(self, fasta_filename):
        with open(fasta_filename, 'r') as infile:
            transcript = infile.readline().strip().split(' ')[0].replace('>', '')
            sequence = ''
            for line in infile:
                if line[0] == '>':
                    self.tran2sequence[transcript] = sequence
                    transcript = line.strip().split(' ')[0].replace('>', '')
                    sequence = ''
                else:
                    sequence += line
            self.tran2sequence[transcript] = sequence
        return
    
    def write_gene_fa(self, storage_prefix, gene):
        output_filename = '{}gene_fastas/{}.fa'.format(storage_prefix, gene.name)
        
        with open(output_filename, 'w') as outfile:
            for tran, seq in gene.trans.items():
                outfile.write('>{} gene:{}\n'.format(tran, gene.name))
                outfile.write('{}\n'.format(seq))
        
        return output_filename
    
    def write_new_fasta(self, full_fasta_filename, fasta_input_filename):
        # create transcript to gene dictionary
        tran2gene = collections.defaultdict(set)
        for gene_name in self.genes:
            for tran in self.genes[gene_name].trans:
                tran2gene[tran].add(gene_name)
        
        # write the full fasta output
        with open(full_fasta_filename, 'w') as outfile:
            for chunk in self.get_fasta_chunks(fasta_input_filename):
                seg = chunk.split('\n', 1)
                transcript = seg[0].strip().split(' ')[0].replace('>', '')
                sequence = seg[1]
                
                outfile.write('>{} gene:{}\n'.format(transcript, ','.join(list(tran2gene[transcript]))))
                outfile.write(sequence + '\n')
        return