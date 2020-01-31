'''
Sean West
30 January 2020
'''
import sys
import getopt

class Input():
    def __init__(self):
        self.args = self.parse_command_line()
        
    def help(self):
        out = '\n'
        out += 'command: transfigg_prep\n'
        out += 'purpose: parse GMAP output to produce distance matrices using MUSCLE, MAFFT, Clustal Omega, and t-coffee\n'
        out += '\n'
        
        out += 'required parameters:\n'
        out += '\t-g  --gmap_output=\tThe exact output file produced by gmap.\n'
        out += '\t-f  --fa=\t\tThe FULL-unmapped de novo transcriptome.\n'
        out += '\t-s  --storage_prefix=\tThe directory that will store all the intermediate and output files.\n'
        out += '\t-c  --cores=\t\tThe number of cores that should be running concurrently when processing the MSA jobs.\n'
        out += '\t-m  --ram=\t\tThe amount of RAM that should be taken up by the concurrent MSA jobs.\n'
        out += '\t-r  --reference=\tThe filename of a reference file that has gene names and genomic locations\n'
        out += '\n'
        
        out += 'help:\n'
        out += '\t-h  --help\tPrint this message.\n'
        out += '\n'
        
        sys.stdout.write(out)
        return
    
    def parse_command_line(self):        
        shortopts = 'g:f:s:c:r:h'
        longopts = ['gmap_output=', 'fa=', 'storage_prefix=', 'cores=', 'ram=', 'help']
        
        try:
            opts, _ = getopt.getopt(sys.argv[1:], shortopts, longopts)
        except getopt.GetoptError as err:
            print(err)
            self.help()
            sys.exit(2)
            
        args = {}
        for o, a in opts:
            if o == '-g' or o == '--gmap_output':
                args['gmap_output'] = a
            elif o == '-f' or o == '--fa':
                args['fa'] = a  
            elif o == '-s' or o == '--storage_prefix':
                if a[-1] != '/':
                    a += '/'
                args['storage_prefix'] = a
            elif o == '-c' or o == '--cores':
                args['cores'] = a
            elif o == '-r' or o == '--ram':
                args['ram'] = a
            elif o == '-h' or o == '--help':
                self.help()
                sys.exit(0)
        
        if not opts:
            self.help()
            sys.exit(0)
            
        return args
        
    