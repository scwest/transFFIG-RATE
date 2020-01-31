'''
This algorithm has room for improvement. 
Perhaps there is a way to check directly after each hierarchical aggregation,
rather than test 100+ points along the hierarchy for differences.
'''

from scipy.cluster import hierarchy

class Clustering():
    def __init__(self, input):
        self.checks = [x/1000 for x in range(1000)]
        self.method = input.args['clustering']
    
    def cluster(self, gene):
        if len(gene.transcripts) <= 1:
            yield gene.transcripts
        elif self.method == 'hierachical':
            for clusters in self.hierarchical(gene):
                yield clusters
            
    def hierarchical(self, gene):
        fit = hierarchy.linkage(gene.distance) # this conducts the whole clustering
        previous = 0
        for check in self.checks:
            clusters = hierarchy.fcluster(fit, check, criterion='distance') # check for new cut at every 1/1000 points
            if previous != len(clusters): # check to make sure there is an actual change in the clustering
                yield clusters
                previous = len(clusters)
