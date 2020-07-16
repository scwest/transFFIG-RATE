'''
This is the main control file. This can be called directly to run the code.
See input.py for a list of inputs, or call this file with "-h"/"--help".
'''

from transffig_rate import Input
from transffig_rate import Transcript
from transffig_rate import Gene
from transffig_rate import Distance
from transffig_rate import Clustering
from transffig_rate import Reliability
from transffig_rate import Trustworthiness
from transffig_rate import Export 
from transffig_rate import Storage
from transffig_rate import Verbose

class Control():
    def main(self):
        input = Input() 
        
        verbose = Verbose(input)
        verbose.welcome()
        
        storage = Storage(input)
        distance = Distance(input)
        clustering = Clustering(input)
        reliability = Reliability()
        trustworthiness = Trustworthiness()
        export = Export(input)
        
        verbose.direct('Uploading Genes')
        genes = storage.upload_genes(Gene, Transcript)
        
        verbose.direct('Running Gene by Gene')
        a = input.args['alpha']
        c = 0
        t = len(genes)
        for gene in genes:
            c += 1
            if len(gene.transcripts) < 3:
                continue
            verbose.mirage('\tGene {} of {} - Uploading Distance Matrix'.format(c, t))
            
            gene.distance = distance.get_matrix(gene)
            
            verbose.mirage('\tGene {} of {} - Setting-up Base Entropies (Trustworthiness)'.format(c, t))
            gene = trustworthiness.get_base_entropy(gene)
            
            verbose.mirage('\tGene {} of {} - Clustering and Calculating Reliabilities/Trustworthiness'.format(c, t))
            try:
                for clusters in clustering.cluster(gene):
                    gene.granularity += [(len(clusters) - 1) / (len(gene.transcripts) - 1)] # granularity as the [0,1] scale of min/max clusters
                    gene.reliability += [reliability.calculate(gene, clusters)]
                    gene.trustworthiness += [trustworthiness.calculate(gene, clusters)]
                    
                    '''
                    print('PRINTING')
                    #print('{} vs {}'.format(gene.rt_max, gene.reliability[-1] + gene.trustworthiness[-1]))
                    print(gene.reliability)
                    print(gene.trustworthiness)
                    print('PRINTING')
                    '''
                    
                    if gene.rt_max < a*gene.reliability[-1] + (1-a)*gene.trustworthiness[-1]:
                        gene.rt_max = a*gene.reliability[-1] + (1-a)*gene.trustworthiness[-1]
                        gene.rt_clusters = clusters
            except:
                print('\tright here <{}>'.format(gene.name))
                raise Exception
            
            verbose.mirage('\tGene {} of {} - Exporting                  '.format(c, t))
            export.gene(gene)
            
            '''
            print('\nSTART')
            print(gene.reliability)
            print(gene.trustworthiness)
            print('END')
            '''
        verbose.mirage_close()
        
        verbose.direct('Full Export')
        verbose.direct('\tGranularities')
        export.set_granularities()
        verbose.direct('\tReliabilities')
        export.set_reliabilities()
        verbose.direct('\tTrustworthiness')
        export.set_trustworthiness()
        verbose.direct('\tFasta Groups')
        export.fastg(storage)
        verbose.direct('\tExpression')
        export.expression(storage)
        
        return
    
def smain():
    stick = Control()
    stick.main()
    return

