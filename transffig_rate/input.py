import sys
import getopt

class Input():
    def __init__(self):
        self.args = self.parse_command_line()
        
    def help(self):
        out = '\n'
        out += 'command: transfigg\n'
        out += 'purpose: organize sub-gene transcripts into functional groups for directed downstream expression analysis\n'
        out += '\n'
        
        out += 'required parameters:\n'
        out += '\t-t  --transcripts=\t\tThe fasta of transcripts (output of de novo assembly)\n'
        out += '\t-g  --mapping=\t\tFile with structure: "<gene>\t<transcript1>,<transcript2>,...,<final transcript>\\n"\n'
        out += '\t-o  --output=\t\tPath and file name of the fasta output\n'
        
        out += 'optional parameters:\n'
        out += '\t-m  --method=\tChoice of method to use for clustering/distance (default: hierarchical/domain_function)\n'
        out += '\t\t\t\tclustering: hierarchical; distance: domain_function,msa\n'
        
        out += 'distance: domain_function'
        out += '\t-d  --domain_prediction\t\tFile with structure: "<transcript>\t<domain1>,<domain2>,...,<final domain>\\n"\n'
        out += '\t-i  --domain_interaction\t\tAn interaction file between domains with structure: "<domain1>\t<domain2>\\n"\n'
        out += '\t-f  --domain_function=\t\tA mapping file between domains and functions in the same format as the previous files\n'
        
        out += 'distance: msa'
        out += '\t-p  --msa_path\t\tPath to a directory of the MSA distance matrices. Each of the matrix must be named "gene_name.txt"\n\n'
        
        out += 'help:'
        out += '\t-h  --help\tPrint this message.\n'
        out += '\n'
        
        sys.stdout.write(out)
        return
    
    def parse_command_line(self):
        sys.stdout.write('Collecting input ... ')
        sys.stdout.flush()
        
        shortopts = 't:g:d:i:f:o:hm:p:'
        longopts = ['transcripts=', 'mapping=', \
                    'domain_prediction=', 'domain_interaction=', 'domain_function',\
                    'output=', 'help', 'method=', '--msa_path=']
        try:
            opts, _ = getopt.getopt(sys.argv[1:], shortopts, longopts)
        except getopt.GetoptError as err:
            print(err)
            self.help()
            sys.exit(2)
            
        args = {'clustering':'hierarchical', 'distance':'water'}
        for o, a in opts:
            if o == '-t' or o == '--transcripts':
                args['transcripts'] = a
            elif o == '-g' or o == '--mapping':
                args['mapping'] = a  
            elif o == '-d' or o == '--domain_prediction':
                args['domain_prediction'] = a
            elif o == '-i' or o == '--domain_interaction':
                args['domain_interaction'] = a
            elif o == '-f' or o == '--domain_function':
                args['domain_function'] = a
            elif o == '-o' or o == '--output':
                args['output'] = a
            elif o == '-m' or o == '--method':
                args['clustering'], args['distance'] = a.split('/')
            elif o == '-p' or o == '--msa_path':
                args['msa_path'] = a
            elif o == '-h' or o == '--help':
                self.help()
                sys.exit(0)
        
        if not opts:
            self.help()
            sys.exit(0)
            
        sys.stdout.write('done\n')
        return args
        
    