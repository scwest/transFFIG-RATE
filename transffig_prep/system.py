'''
Sean West
30 January 2020

This portion of the code handles calls to the operating system
to make sure that runs are within the proper core/RAM constraints.
'''

import sys 

class System():
    def __init__(self, ram, cores):
        self.ram = ram
        self.cores = cores
        
        