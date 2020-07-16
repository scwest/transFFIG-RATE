from skbio.alignment import StripedSmithWaterman as SSW
import itertools
import collections
import os
import sys

class Distance():
    def __init__(self, input):
        self.distance_measure = input.args['distance']
        self.matrices_path = input.args['msa_path']
    
    def get_matrix(self, gene):
        if self.distance_measure == 'water':
            return self.water(gene)
        elif self.distance_measure == 'msa':
            return self.msa(gene)
        return []
    
    def msa(self, gene):
        distance_matrix = collections.defaultdict(dict)
        max_distance = 0
        
        distance_matrix_filename = self.matrices_path+gene.name+'.muscle.csv'
        if not os.path.isfile(distance_matrix_filename):
            print('\n\nThe filename was not found:\t{}'.format(distance_matrix_filename))
            print('Gene:\t{}'.format(gene.name))
            #print('Transcripts:\t{}'.format(','.join(list(gene.transcript_dictionary.keys()))))
            print('len(gene.transcripts):\t{}'.format(len(gene.transcripts)))
            print('\n')
            sys.exit(2)
        
        with open(distance_matrix_filename, 'r') as dfile:
            header = dfile.readline().strip().split(',')
            
            for line in dfile:
                line = line.strip().split(',')
                tran1_name = line[0]
                values = [float(x) for x in line[1:]]
                
                for i in range(len(values)):
                    a, b = sorted([tran1_name, header[i]])
                    distance_matrix[a][b] = values[i]
        
        dm = []
        transcripts = sorted(list(distance_matrix.keys()))
        for transcript1 in transcripts[:]:
            temp = []
            for transcript2 in transcripts[:]:
                a, b = sorted([transcript1, transcript2])
                temp.append(distance_matrix[a][b])
            dm.append(temp)
                    
        return dm
    
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