


class Gene():
    def __init__(self, name):
        self.transcripts = []
        self.distance = []
        self.reliability = []
        self.granularity = []
        self.trustworthiness = []
        self.name = name
        self.expression = []
        self.rt_max = 0
        self.rt_clusters = []
        self.transcript_dictionary = {}
        
    def calculate_expression(self):
        self.expression = self.transcripts[0].expression
        for transcript in self.transcripts[1:]:
            self.expression = [x+y for x,y in zip(self.expression, transcript.expression)]
        for i in range(len(self.expression)):
            if not self.expression[i]:
                self.expression[i] = 1
        return