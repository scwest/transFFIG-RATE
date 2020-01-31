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


class Control():
    def main(self):
        input = Input() 
        
        storage = Storage(input)
        distance = Distance(input)
        clustering = Clustering()
        reliability = Reliability()
        trustworthiness = Trustworthiness()
        export = Export(input)
        
        genes = storage.upload(Gene, Transcript)
        for gene in genes:
            gene.distance = distance.get_matrix(gene)
            gene = trustworthiness.get_base_entropy(gene)
            
            for clusters in clustering.cluster(gene, storage):
                gene.granularity += [(len(clusters) - 1) / (len(gene.transcripts) - 1)] # granularity as the [0,1] scale of min/max clusters
                gene.reliability += [reliability.calculate(gene, clusters)]
                gene.trustworthiness += [trustworthiness.calculate(gene, clusters)]
                
                if gene.rt_max < gene.reliability[-1] + gene.trustworthiness[-1]:
                    gene.rt_clusters = clusters
                    
            export.gene(gene)
        
        outdir = input.args['outdir']
        prefix = input.args['prefix']
        export.granularities(outdir+prefix+'granularities.txt')
        export.reliabilities(outdir+prefix+'reliabilities.txt')
        export.trustworthiness(outdir+prefix+'trustworthiness.txt')
        export.fastg(outdir+prefix+'transcript_groups.fastg')
        export.expression(outdir+prefix+'expression.csv')
        
        return
    
def smain():
    stick = Control()
    stick.main()
    return

stick = Control()
stick.main()
