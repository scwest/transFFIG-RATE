'''
This algorithm has room for improvement. 
Perhaps there is a way to check directly after each hierarchical aggregation,
rather than test 100+ points along the hierarchy for differences.
'''

from scipy.cluster import hierarchy
import collections

class Clustering():
    def __init__(self, input):
        self.checks = [x/1000 for x in range(1000)]
        self.method = input.args['clustering']
        
    
    def _to_names(self, clusters, transcripts):
        '''
        print('\nhere')
        print(clusters)
        print(transcripts)
        print('here\n')
        '''
        dclusters = collections.defaultdict(list)
        for i in range(len(clusters)):
            dclusters[clusters[i]-1] += [transcripts[i]]
            
        nclusters = []
        for k in sorted(list(dclusters.keys())):
            nclusters.append(dclusters[k])
        
        return nclusters
    
    def cluster(self, gene):
        if len(gene.transcripts) <= 1:
            yield gene.transcripts
        elif self.method == 'hierarchical':
            for c in self.hierarchical(gene):
                yield c
            '''
            #This will forward the generator anyway:
            for clusters in self.hierarchical(gene):
                yield clusters
            '''
        return
            
    def hierarchical(self, gene):
        fit = hierarchy.linkage(gene.distance) # this conducts the whole clustering
        previous = len(gene.transcripts)
        yield self._to_names(list(range(1, len(gene.transcripts)+1)), sorted([x.name for x in gene.transcripts]))
        for check in self.checks:
            clusters = hierarchy.fcluster(fit, check, criterion='distance') # check for new cut at every 1/1000 points
            clusters = self._to_names(clusters, sorted([x.name for x in gene.transcripts]))
            if previous != len(clusters): # check to make sure there is an actual change in the clustering
                yield clusters
                previous = len(clusters)
        if previous != 1:
            yield self._to_names([1] * len(gene.transcripts), sorted([x.name for x in gene.transcripts]))
                
                
