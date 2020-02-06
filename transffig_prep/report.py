'''
Sean West
5 February 2020
'''
import collections
import numpy as np

class Report():
    def run(self, genes):
        print('Number of Genes:\t{}'.format(len(genes)))
        
        c = sum([1 for gene in genes if len(genes[gene].trans) > 2])
        print('Number of Genes with > 2 Transcripts:\t{}'.format(c))
        
        tran2genes = collections.defaultdict(set)
        for gene_name in genes:
            for tran in genes[gene_name].trans:
                tran2genes[tran].add(gene_name)
        print('Number of Transcripts:\t{}'.format(len(tran2genes)))
        print('Number of Transcripts with > 1 gene:\t{}'.format(len([x for x in tran2genes if len(tran2genes[x]) > 1])))
        
        glengths = []
        for gene_name in genes:
            glength = abs(genes[gene_name].start - genes[gene_name].end)
            glengths.append(glength)
        print('Mean Gene Length:\t{}'.format(np.mean(glengths)))
        
        tlengths = []
        for tran_name in tran2genes:
            tlength = len(genes[list(tran2genes[tran_name])[0]].trans[tran_name])
            tlengths.append(tlength)
        print('Mean Transcript Length:\t{}'.format(np.mean(tlengths)))
        print('Min Transcript Length:\t{}'.format(min(tlengths)))
        
        return
