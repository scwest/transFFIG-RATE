'''
Sean West
5 February 2020

GMAP output is separated into a long text list of queries and
their associated hits. Each of these chunks corresponds to
a sequence query using in GMAP. This module parses out this text
into meaningful values, mapping transcripts to genomic locations.
'''

class Gmap_chunk():
    def __init__(self, text):
        self.text = text
        self.hits = []
        
    class Hit():
        def __init__(self):
            self.name = ''
            self.chromosome = ''
            self.start = 0
            self.end = 0
            self.strand = 1
        
    def get_transcript(self):
        return self.text.split('\n')[0].split(' ')[0].replace('>', '')
    
    def get_sequence(self):
        return '\n'.join(self.text.strip().split('\n')[1:])
        
    def process(self):
        inpaths = False
        fid = self.text.split('\n', 1)[0].split(' ')[0].replace('>', '')
        for line in self.text.split('\n'):
            if line[:5] == 'Paths':
                inpaths = True
                continue
            elif line[:10] == 'Alignments':
                inpaths = False
                continue
            if inpaths:
                if line[:6] == '  Path':
                    try:
                        loc = line.split(' ')[-3]
                        loc = loc.split(':')
                        chromosome = loc[0]
                        spots = loc[1].split('..')
                        start = int(spots[0].replace(',', ''))
                        end = int(spots[1].replace(',', ''))
                        strand = 1
                        if line.split('(')[-1][0] == '-':
                            strand = -1
                    except:
                        print(line)
                        print(loc)
                        raise
            
                    hit = self.Hit()
                    hit.name = fid
                    hit.chromosome = chromosome
                    hit.start = int(start)
                    hit.end = int(end)
                    hit.strand = strand
                    self.hits.append(hit)
        return
                