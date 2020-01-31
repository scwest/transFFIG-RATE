from skbio.alignment import StripedSmithWaterman as SSW
import itertools
import collections

class Distance():
    def __init__(self, input):
        self.distance_measure = input.args['distance']
    
    def get_matrix(self, gene):
        if self.distance_measure == 'water':
            return self.water(gene)
        elif self.distance_measure == '<placeholder>':
            return []
        return []
    
    def water(self, gene):
        distance_matrix = collections.defaultdict(dict)
        max_distance = 0
        
        for tran1, tran2 in itertools.combinations(gene.transcripts, 2):
            query = SSW(tran1.sequence, score_only=True)
            alignment = query(tran2.sequence)
            distance = alignment['optimal_alignment_score']
            distance_matrix[tran1.name][tran2.name] = distance
            
            if distance > max_distance:
                max_distance = distance
                
        for tran1, tran2 in itertools.combinations(gene.transcripts, 2):
            distance_matrix[tran1.name][tran2.name] = 1 - (distance_matrix[tran1.name][tran2.name] / max_distance)
        return distance_matrix