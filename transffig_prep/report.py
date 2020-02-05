'''
Sean West
5 February 2020
'''
import collections

class Report():
    def run(self, genes, commands):
        print('Number of Genes:\t{}'.format(len(genes)))
        
        c = sum([1 for gene in genes if len(gene.trans) > 2])
        print('Number of Genes with > 2 Transcripts:\t{}'.format(c))
        
        tran2genes = collections.defaultdict(set)
        for gene in genes:
            for tran in gene.trans:
                tran2genes[tran].add(gene)
        print('Number of Transcripts:\t{}'.format(len(tran2genes)))
        print('Number of Transcripts with > 1 gene:\t{}'.format(len([x for x in tran2genes if len(tran2genes[x]) > 1]))) 
        
        return
