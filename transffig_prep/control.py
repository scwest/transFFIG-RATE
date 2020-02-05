'''
Sean West
30 January 2020

'''

from transffig_prep import Input
from transffig_prep import System
from transffig_prep import Gmap
from transffig_prep import Msa


class Control():
    def main(self):
        print('Reading Inputs')
        input = Input()
        inputs = input.args
        
        # process gmap to create a list of commands that need to be written
        # while putting all the tables for each gene in <storage_prefix/gene_fastas/{gene}.fa>
        print('Checking for existing parsed GMAP output / commands.')
        gmap = Gmap()
        print('Parsing GMAP output (this will take a while).')
        commands = gmap.parse(commands, inputs['storage_prefix'], inputs['fasta'], inputs['gmap_output'], inputs['reference'])
        
        # FOR TESTING ONLY
        
        
        # setup system constraints for running the msa jobs
        print('Setting up system constraints.')
        system = System(inputs['ram'], inputs['cores'])
        
        # run the msa jobs
        # this will check the system to make sure it can continue to run them.
        msa = Msa(system)
        #print('Checking for previously run commands.')
        #commands = msa.check_previous(commands, inputs['storage_prefix'])
        print('Running remaining commands.')
        msa.run_all_commands(commands)
        
        print('Finished. You will find any distance matrices under your listed storage_prefix/distance_matrices/<msa_method>/<gene>.csv.')
        
        return
    
def smain():
    stick = Control()
    stick.main()
    return

