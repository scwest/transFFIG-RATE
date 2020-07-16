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
        vars = {}
        for transcript in gene.transcripts:
            vars[transcript.name] = np.var(transcript.expression)
            
        sum_var = sum(vars.values())
        if sum_var == 0:
            return 1
        
        mrv = []
        perfects = 0
        for x in vars.values():
            if x == 0:
                perfects += 1
            else:
                mrv.append(1/x)
        if perfects:
            mrv += [np.max(mrv)]*perfects
        mean_rec_var = np.mean(mrv)
        
        rec_sum_var = 1/sum_var
        
        per_cluster = 0
        pcs = []
        perfects = 0
        for cluster in clusters:
            a = sum([vars[transcript_name] for transcript_name in cluster])
            if a == 0:
                perfects += 1
            else:
                pcs.append(1/a)
        if perfects:
            pcs += [np.max(pcs)]*perfects
        per_cluster = np.mean(pcs)
        
        if mean_rec_var < rec_sum_var:
            if mean_rec_var - rec_sum_var < 10**-10:
                return mean_rec_var
            else:
                print('Why is the max less than the min!!!:???!?')
                print(list(vars.values()))
                print(mean_rec_var)
                print(rec_sum_var)
                raise Exception
        
        
        return (mean_rec_var - per_cluster) / (mean_rec_var - rec_sum_var)