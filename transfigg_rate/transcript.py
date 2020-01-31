from skbio import DNA

class Transcript():
    def __init__(self):
        self.sequence = DNA('')
        self.name = ''
        #self.expression = [] # I don't think we will need the expression for this code.