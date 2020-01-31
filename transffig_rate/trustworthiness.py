'''
This equation is from the extended PHD Proposal of Sean West.
This section can be found in docs/.
'''
import numpy as np

class Trustworthiness():
    def __init__(self):
        pass
    
    def calculate(self, gene, clusters):
        p = max(clusters)
        
        expressions = []
        for cnum in range(1, p+1):
            expressions.append(self.aggregate([gene.transcripts[x].expression for x in range(gene.transcripts) if clusters[x] == cnum]))
        
        entropies = [self.entropy(x, gene.expression)/np.log2(p) for x in expressions]
        return 1 - np.abs(np.mean(entropies) - gene.base_mean_entropy)
    
    def aggregate(self, vectors):
        aggregation = vectors.pop()
        for vector in vectors:
            aggregation = [x+y for x,y in zip(aggregation, vector)]
        return aggregation
    
    def entropy(self, expression, total_expression):
        ratios = [expression[i]/total_expression[i] for i in range(expression)]
        return -sum([x*np.log2(x) for x in ratios])
    
    def get_base_entropy(self, gene):
        max_entropy = np.log2(len(gene.transcripts))
        entropies = [self.entropy(gene.transcripts[i].expression, gene.expression)/max_entropy for i in range(gene.transcripts)]
        gene.base_mean_entropy = np.mean(entropies)
        return gene