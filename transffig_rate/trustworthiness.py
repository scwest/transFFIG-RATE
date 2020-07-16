'''
This equation is from the extended PHD Proposal of Sean West.
This section can be found in docs/.
'''
import numpy as np

class Trustworthiness():
    def __init__(self):
        pass
    
    def calculate(self, gene, clusters):
        if len(clusters) == 1:
            return 0
        
        entropies = []
        max_entropy = np.log2(len(clusters))
        for sample_number in range(len(gene.transcripts[0].expression)):
            per_cluster_expression = []
            for cluster in clusters:
                expressions = [gene.transcript_dictionary[transcript_name].expression[sample_number] for transcript_name in cluster]
                per_cluster_expression.append(sum(expressions))
            sper = sum(per_cluster_expression)
            if sper == 0:
                return 1
            entropy = self.entropy(per_cluster_expression, sper) / max_entropy
            entropies.append(entropy)
        
        '''
        for cluster in clusters:    
            cluster_expressions = [gene.transcript_dictionary[transcript].expression for transcript in cluster]
            agg = self.aggregate(cluster_expressions)
            expressions.append(agg)
        
        entropies = [self.entropy(x, gene.expression)/np.log2(len(clusters)) for x in expressions]
        
        for entropy in entropies:
            if entropy < 0 or entropy > 1:
                print('BASE ENTROPY OUT OF 0,1')
                print(entropies)
                print('')
                raise Exception
        '''
        
        return 1 - np.abs(np.mean(entropies) - gene.base_mean_entropy)
    
    def aggregate(self, vectors):
        aggregation = vectors.pop()
        for vector in vectors:
            aggregation = [x+y for x,y in zip(aggregation, vector)]
        return aggregation
    
    def entropy(self, expression, total_expression):
        ratios = [expression[i]/total_expression for i in range(len(expression))]
        entropy = -sum([x*np.log2(x) for x in ratios if x])
        return entropy
    
    def get_base_entropy(self, gene):
        max_entropy = np.log2(len(gene.transcripts))
        
        # by sample
        entropies = []
        for sample_number in range(len(gene.transcripts[0].expression)):
            expression = [transcript.expression[sample_number] for transcript in gene.transcripts]
            sexpr = sum(expression)
            if sexpr == 0:
                gene.base_mean_entropy = 1
                return gene
            entropy = self.entropy(expression, sexpr) / max_entropy
            entropies.append(entropy)
            
            if entropy < 0 or entropy > 1:
                print('BASE ENTROPY OUT OF 0,1')
                print(max_entropy)
                print(entropies)
                print([x*max_entropy for x in entropies])
                print(gene.transcripts[0].expression)
                print(gene.expression)
                print('')
                raise Exception
            
        gene.base_mean_entropy = np.mean(entropies)
        return gene