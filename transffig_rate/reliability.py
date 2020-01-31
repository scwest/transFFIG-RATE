'''
This equation is from the extended PHD Proposal of Sean West.
This section can be found in docs/.
This is not true statistical reliability, rather it is the ratio of:
  (observed reliability - min reliability) / (max reliability - min reliability)
  which removes the variance of the error from the equation.
  Estimation of the variance of the error is time consuming and unnecessary for 
  finding relative reliabilities, IF!!! we can assume that the variance of the error
  is equal across all the transcripts, which we do. See the paper in docs/ for more details.
'''

import numpy as np

class Reliability():
    def __init__(self):
        pass
    
    def calculate(self, gene, clusters):
        vars = [np.var(x.expression) for x in gene.transcripts]
        mean_rec_var = np.mean([1/x for x in vars])
        rec_sum_var = 1/sum(vars)
        per_cluster = 0
        max_cnum = max(clusters)
        for i in range(1, max_cnum+1):
            per_cluster += 1/sum([vars[x] for x in range(len(clusters)) if clusters[x] == i])
        per_cluster = per_cluster / max_cnum
        
        return (mean_rec_var - per_cluster) / (mean_rec_var - rec_sum_var)