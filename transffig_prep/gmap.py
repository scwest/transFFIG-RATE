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

class Gmap():
    def __init__(self):
        self.genes = {}
        self.current_gene_number = 1
        
    class Gene():
        def __init__(self, chromosome='', start=0, end=0, trans=collections.defaultdict(str), name='', strand=1):
            self.chromosome = chromosome
            self.start = start
            self.end = end
            self.trans = trans
            self.fa_filename = ''
            self.name = name
            self.strand = strand
            
    class Chunk():
        def __init__(self, text):
            self.text = text
            
        def get_transcript(self):
            return self.text.split('\n')[0].split(' ')[0].replace('>', '')
        
        def get_sequence(self):
            return '\n'.join(self.text.strip().split('\n')[1:])
            
        def process(self):
            inpaths = False
            fid = self.text.split('\n', 1)[0].split(' ')[0].replace('>', '')
            for line in self.text.split('\n'):
                if line[:5] == 'Paths':
                    inpaths = True
                    continue
                elif line[:10] == 'Alignments':
                    inpaths = False
                    continue
                if inpaths:
                    if line[:6] == '  Path':
                        try:
                            loc = line.split(' ')[-3]
                            loc = loc.split(':')
                            chromosome = loc[0]
                            spots = loc[1].split('..')
                            start = int(spots[0].replace(',', ''))
                            end = int(spots[1].replace(',', ''))
                            strand = 1
                            if line.split('(')[-1][0] == '-':
                                strand = -1
                        except:
                            print(line)
                            print(loc)
                            raise
                        yield fid, chromosome, start, end, strand
            return
                        
    def absorb_chunk(self, tran, chromosome, start, end, strand):
        '''
        New genes receive an 'arbitrary_<>' gene id.
        This is to emphasize that they are not real ids. That they may be reassigned to real
          gene ids if desired.
          
        ABSORB CHUNK DOES NOT ASSIGN SEQUENCES TO TRANSCRIPTS
        '''
        overlapped_gids = []
        
        # Check to see all the genes that overlap the genomic location of the transcript
        for gid, gene in self.genes.items():
            if gene.chromosome != chromosome:
                continue
            if gene.strand != strand:
                continue
            if self.isoverlap(gene, start, end):
                overlapped_gids.append(gid)
            
        # Incorporate the transcript into existing (or new) genes
        ### If there is no gene that overlaps, we'll create a new gene
        ### create a new gene at this location
        if len(overlapped_gids) == 0:
            t = collections.defaultdict(str)
            t[tran] = '' # NOTE: that if the gene is constructed this way, the transcripts don't have sequences
            self.genes['arbitrary_{}'.format(self.current_gene_number)] = self.Gene(chromosome=chromosome,\
                                                                                    start=start,\
                                                                                    end=end,\
                                                                                    trans = t,
                                                                                    strand = strand)
            self.current_gene_number += 1
        
        ### If we found 1 gene that the transcript overlaps;
        ### add to an existing gene
        elif len(overlapped_gids) == 1:
            gid = overlapped_gids.pop()
            self.genes[gid].start = min([self.genes[gid].start, start])
            self.genes[gid].end = max([self.genes[gid].end, end])
            self.genes[gid].trans[tran] = ''
        
        ### If we found multiple genes that the transcript overlaps;
        ### These are either genes that are on the same genomic location
        ### or they are the same gene. Either way, add them together.
        else:
            gene = self.Gene()
            gene.start = min([self.genes[x].start for x in overlapped_gids]+[start])
            gene.end = max([self.genes[x].end for x in overlapped_gids]+[end])
            gene.chromosome = self.genes[overlapped_gids[0]].chromosome # they're all the same; only [0],[1] are guaranteed
            gene.strand = self.genes[overlapped_gids[0]].strand
            for g in [self.genes[x] for x in overlapped_gids]:
                for t, seq in g.trans.items():
                    gene.trans[t] = seq
            
            self.genes['arbitrary_{}'.format(self.current_gene_number)] = gene
            self.current_gene_number += 1
            
            for i in overlapped_gids:
                del self.genes[i]
        return
        
    def check(self, storage_prefix):
        commands = []
        storage_filename = '{}.commands.txt'.format(storage_prefix)
        
        if os.path.isfile(storage_filename):
            with open(storage_filename, 'r') as infile:
                for line in infile:
                    commands.append(line.strip().split(' '))
        return commands
    
    def command(self, gene, storage_prefix):
        commands = []
        for call in ['transffig_muscle', 'transffig_mafft', 'transffig_tcoffee', 'transffig_clustalo']:
            commands.append([call, gene.fa_filename, '{}distance_matrices/{}.csv'.format(gene.name, storage_prefix)])
        
        return commands
    
    def get_chunks(self, input_filename):
        with open(input_filename, 'r') as infile:
            chunk = infile.readline()
            for line in infile:
                if line[0] == '>':
                    yield self.Chunk(chunk) 
                    chunk = line
                else:
                    chunk += line
        return self.Chunk(chunk)
    
    def isoverlap(self, gene, start, end):
        if start > gene.start and start < gene.end:
            return True
        if end > gene.start and end < gene.end:
            return True
        return False
    
    def link_genes(self, reference_filename):
        '''
        The file must be setup as:
        <gene name>,<chromosome>,<start>,<end>
        '''
        ref = []
        with open(reference_filename, 'r') as infile:
            for line in infile:
                ref = line.strip().split(',')
                ref[2] = int(ref[2])
                ref[3] = int(ref[3])
                ref[4] = int(ref[4])
        
        for gene_name, gene in self.genes.items():
            for rgene in ref:
                if ref[1] != gene.chromosome:
                    continue
                if ref[4] != gene.strand:
                    continue
                if isoverlap(gene, rgene[1], rgene[2]):
                    self.genes[ref[0]] = gene
                    self.genes[ref[0]].name = ref[0]
                    del self.genes[gene_name]
                    break
        return
    
    def parse(self, commands, storage_prefix, fasta_input_filename, gmap_output_filename, reference_filename):
        commands = []
        fasta_output_filename = '{}full.fa'.format(storage_prefix)
        
        # setup output directories
        fasta_dir = '{}gene_fastas'.format(storage_prefix)
        subprocess.call(['mkdir', fasta_dir])
        dm_dir = '{}distance_matrices'.format(storage_prefix)
        subprocess.call(['mkdir', dm_dir])
        
        # check for an existing .fa that has the gmap already processed
        if os.path.isfile(fasta_output_filename):
            print('\tParsing Existing Commands.')
            self.parse_full_fa(fasta_output_filename, storage_prefix)
        else:
            print('\tParsing GMAP hits (chunks)')
            self.parse_gmap(gmap_output_filename, fasta_input_filename, fasta_output_filename)
            
            # if there is a reference fasta, link the genes to their real names (referenced names)
            print('\tLinking genes to the given reference.')
            if reference_filename:
                self.link_genes(reference_filename)
                
            # write the genes to their individual fa's
            print('\tWriting gene-specific fasta files.')
            for gene_name, gene in self.genes.items():
                self.genes[gene_name].fa_filename = self.write_gene_fa(storage_prefix, gene)
        
        # write the commands for each gene
        for gene in self.genes.values():
            commands += self.command(gene, storage_prefix)
        
        return commands 
    
    def parse_full_fa(self, fasta_output_filename, storage_prefix):
        with open(fasta_output_filename, 'r') as infile:
            for line in infile:
                if line[0] == '>':
                    tran_name = line.strip().split(' ')[0].replace('>', '')
                    gene_name = line.strip().split('gene:')[-1].split(' ')[0]
                    
                    if gene_name not in self.genes:
                        self.genes[gene_name] = self.Gene()
                        self.genes[gene_name].name = gene_name
                        self.genes[gene_name].fa_filename = '{}/gene_fastas/{}.fa'.format(storage_prefix, gene_name)
                else:
                    self.genes[gene_name].trans[tran_name] += line
        return
    
    def parse_gmap(self, gmap_output_filename, fasta_input_filename, fasta_output_filename):
        # get all the genes
        c = 1
        for chunk in self.get_chunks(gmap_output_filename):
            sys.stdout.write('\r\t\tChunk: {}'.format(c))
            sys.stdout.flush()
            c += 1
            for tran, chromosome, start, end, strand in chunk.process():
                self.absorb_chunk(tran, chromosome, start, end, strand)
        print('')
        
        # create transcript to gene dictionary
        tran2gene = collections.defaultdict(set)
        for gene_name in self.genes:
            for tran in self.genes[gene_name].trans:
                tran2gene[tran].add(gene_name)
        
        # write the full fasta output
        multiples = 0 # the number of transcripts that are mapped to multiple genes
        with open(fasta_output_filename, 'w') as outfile:
            for chunk in self.get_chunks(fasta_input_filename):
                transcript = chunk.get_transcript()
                sequence = chunk.get_sequence()
                gene_names = tran2gene[transcript]
                
                if len(gene_names) > 1:
                    multiples += 1
                if not gene_names:
                    continue
                    
                for gene_name in gene_names:
                    if gene_name in self.genes:
                        self.genes[gene_name].trans[transcript] = sequence
                        break
                
                outfile.write('>{} gene:{}\n'.format(transcript, gene_name))
                outfile.write(sequence + '\n')
                    
        return
    
    def write_gene_fa(self, storage_prefix, gene):
        output_filename = '{}gene_fastas/{}.fa'.format(storage_prefix, gene.name)
        
        with open(output_filename, 'w') as outfile:
            for tran, seq in gene.trans.items():
                outfile.write('>{} gene:{}\n'.format(tran, gene.name))
                outfile.write('{}\n'.format(seq))
        
        return output_filename
