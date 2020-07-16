import collections

class Storage():
    def __init__(self, input):
        self.args = input.args
        self.sequences = {}
        self.sample_ids = []
    
    def upload_genes(self, Gene, Transcript):
        genes = []
        
        transcript2expression = self._upload_expression()
        gene2transcripts = self._upload_mapping()
        
        # add code to upload the mapping file and the expression files
        for gene_name, transcript_names in gene2transcripts.items():
            if gene_name == '':
                print('How did this happen?')
                raise Exception
                continue
            gene = Gene(gene_name)
            for transcript_name in transcript_names:
                transcript = Transcript(transcript_name)
                transcript.expression = transcript2expression[transcript_name]
                gene.transcripts.append(transcript)
                gene.transcript_dictionary[transcript_name] = transcript
            gene.calculate_expression()
            genes.append(gene)
        
        return genes

    def _upload_expression(self):
        transcript2expression = collections.defaultdict(list)
        
        with open(self.args['expression'], 'r') as efile:
            header = efile.readline()
            self.sample_ids = header.strip().split(',')[1:]
            for line in efile:
                line = line.strip().split(',')
                transcript = line[0]
                expression = [float(x) for x in line[1:]]
                transcript2expression[transcript] = expression
                
        return transcript2expression
    
    def _upload_mapping(self):
        gene2transcripts = collections.defaultdict(list)
        
        if 'mapping' in self.args:
            with open(self.args['mapping'], 'r') as mfile:
                for line in mfile:
                    gene_name, transcript_names = line.strip().split('\t')
                    transcript_names = transcript_names.split(',')
                    gene2transcripts[gene_name] = list(set(transcript_names))
        else:
            with open(self.args['transcripts'], 'r') as infile:
                line = infile.readline()
                line = line.strip().split(' ')
                
                sequence = ''
                gene_names = line[1].replace('gene:', '').split(',')
                transcript = line[0].replace('>', '')
                for gene_name in gene_names:
                    gene2transcripts[gene_name] += [transcript]
                 
                for line in infile:
                    if line[0] == '>':
                        if sequence and transcript not in self.sequences and gene_names[0] != '':
                            self.sequences[transcript] = sequence
                        sequence = ''
                        
                        line = line.strip().split(' ')
                        gene_names = line[1].replace('gene:', '').split(',')
                        transcript = line[0].replace('>', '')
                        for gene_name in gene_names:
                            gene2transcripts[gene_name] += [transcript]
                        
                    else:
                        sequence += line.strip()
                
                if sequence and transcript not in self.sequences and gene_names[0] != '':
                    self.sequences[transcript] = sequence
                
        return gene2transcripts
    
    