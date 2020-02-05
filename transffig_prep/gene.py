'''
Sean West
5 February 2020
'''

class Gene():
    def __init__(self, chromosome='', start=0, end=0, trans=collections.defaultdict(str), name='', strand=1):
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.trans = trans
        self.fa_filename = ''
        self.name = name
        self.strand = strand
        
    