import collections


class Export():
    def __init__(self):
        self.ingrans = collections.defaultdict(set)
        self.ingran_expression = collections.defaultdict(list)
        self.reliabilities = collections.defaultdict(list)
        self.trustworthiness = collections.defaultdict(list)
        self.granularities = collections.defaultdict(list)
        self.all_genes = []
    
    def gene(self, gene):
        self.granularities[gene.name] = gene.granularities
        self.trustworthiness[gene.name] = gene.trustworthiness
        self.reliabilities[gene.name] = gene.reliabilies
        self.get_ingran(gene)
        self.all_genes.append(gene.name)        
        return
    
    def get_ingran(self, gene):
        max_cnum = max(gene.clusters)
        for i in range(1, max_cnum+1):
            name = 'gene_'+gene.name+'_trangroup_'+i
            
            tran_expressions = [gene.transcripts[x].expression for x in range(len(gene.transcripts)) if gene.clusters[x] == i]
            tran_names = [gene.transcripts[x].names for x in range(len(gene.transcripts)) if gene.clusters[x] == i]
            
            expression = tran_expression.pop()
            for tex in tran_expressions:
                expression = [x+y for x,y in zip(expression, tex)]
                
            self.igran_expression[name] = expression
            self.igrans[name] = set(tran_names)
        return
    
    def granularities(self, outputfilename):
        gene_names = sorted(self.all_genes)
        with open(outputfilename, 'w') as outfile:
            for gene_name in gene_names:
                outfile.write(gene_name + '\t' + ','.join([str(x) for x in self.granularities[gene_name]]) + '\n')
        return
    
    def reliabilities(self, outputfilename):
        gene_names = sorted(self.all_genes)
        with open(outputfilename, 'w') as outfile:
            for gene_name in gene_names:
                outfile.write(gene_name + '\t' + ','.join([str(x) for x in self.reliabilities[gene_name]]) + '\n')
        return
    
    def trustworthiness(self, outputfilename):
        gene_names = sorted(self.all_genes)
        with open(outputfilename, 'w') as outfile:
            for gene_name in gene_names:
                outfile.write(gene_name + '\t' + ','.join([str(x) for x in self.trustworthiness[gene_name]]) + '\n')
        return
    
    def fastg(self, storage, outputfilename):
        ingran_names = sorted(list(self.ingrans.keys()))
        with open(outputfilename, 'w') as outfile:
            for ingran_name in ingran_names:
                outfile.write('>'+ingran_name+'\n')
                for tran in self.ingrans[ingran_name]:
                    outfile.write('>>'+tran+'\n')
                    outfile.write(storage.sequences[tran] + '\n')
        return
    
    def expression(self, storage, outputfilename):
        ingran_names = sorted(list(self.ingran_expression.keys()))
        with open(outputfilename, 'w') as outfile:
            outfile.write(',' + ','.join(storage.sample_ids) + '\n')
            for igran_name in igran_names:
                outfile.write(igran_name + ',' + ','.join([str(x) for x in self.igran_expression[igran_name]]) + '\n')
        return
    
    