import collections
import subprocess
import os

class Export():
    def __init__(self, input):
        self.ingrans = collections.defaultdict(set)
        self.ingran_expression = collections.defaultdict(list)
        self.reliabilities = collections.defaultdict(list)
        self.trustworthiness = collections.defaultdict(list)
        self.granularities = collections.defaultdict(list)
        self.all_genes = []
        self.args = input.args
        self._set_output_prefix()
    
    def _set_output_prefix(self):
        if not os.path.isdir(self.args['output']):
            subprocess.call(['mkdir', self.args['output']])
        return
    
    def gene(self, gene):
        self.granularities[gene.name] = gene.granularity
        self.trustworthiness[gene.name] = gene.trustworthiness
        self.reliabilities[gene.name] = gene.reliability
        self.get_ingran(gene)
        self.all_genes.append(gene.name)        
        return
    
    def get_ingran(self, gene):
        spot = 1
        for cluster in gene.rt_clusters:
            #tran_names = '-'.join([x for x in cluster])
            cluster_name = 'gene_{}_ffig_{}'.format(gene.name, spot)
            spot += 1
            
            tran_expressions = [gene.transcript_dictionary[transcript].expression for transcript in cluster]
            expression = tran_expressions.pop()
            for tex in tran_expressions:
                expression = [x+y for x,y in zip(expression, tex)]
                
            self.ingran_expression[cluster_name] = expression
            self.ingrans[cluster_name] = set(cluster)
        return
    
    def set_granularities(self):
        gene_names = sorted(self.all_genes)
        with open(self.args['output']+'granularities.txt', 'w') as outfile:
            for gene_name in gene_names:
                outfile.write(gene_name + '\t' + ','.join([str(x) for x in self.granularities[gene_name]]) + '\n')
        return
    
    def set_reliabilities(self):
        gene_names = sorted(self.all_genes)
        with open(self.args['output']+'reliabilities.txt', 'w') as outfile:
            for gene_name in gene_names:
                outfile.write(gene_name + '\t' + ','.join([str(x) for x in self.reliabilities[gene_name]]) + '\n')
        return
    
    def set_trustworthiness(self):
        gene_names = sorted(self.all_genes)
        with open(self.args['output']+'trustworthiness.txt', 'w') as outfile:
            for gene_name in gene_names:
                outfile.write(gene_name + '\t' + ','.join([str(x) for x in self.trustworthiness[gene_name]]) + '\n')
        return
    
    def fastg(self, storage):
        ingran_names = sorted(list(self.ingrans.keys()))
        with open(self.args['output']+'transcript_groups.fastg', 'w') as outfile:
            for ingran_name in ingran_names:
                outfile.write('>'+ingran_name+'\n')
                for tran in self.ingrans[ingran_name]:
                    outfile.write('>>'+tran+'\n')
                    if tran in storage.sequences:
                        outfile.write(storage.sequences[tran] + '\n')
        return
    
    def expression(self, storage):
        ingran_names = sorted(list(self.ingran_expression.keys()))
        with open(self.args['output']+'expression.csv', 'w') as outfile:
            outfile.write(',' + ','.join(storage.sample_ids) + '\n')
            for ingran_name in ingran_names:
                outfile.write(ingran_name + ',' + ','.join(['{0:.6f}'.format(x) for x in self.ingran_expression[ingran_name]]) + '\n')
        return
    
    