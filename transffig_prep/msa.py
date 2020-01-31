'''
Sean West
30 January 2020

This code runs the MSA calls within the given system constraints. 
Right now, it only works with the number of cores passed. :(
Implementation of RAM requirements are predictive and difficult to see coming.
I'm not sure if we should ever implement them, but it's possible if we're
  having difficulties.

'''

import subprocess
import os
import time


class Msa():
    def __init__(self, system):
        self.system = system
        
    def check_previous(self, commands, storage_prefix):
        commands_filename = '{}commands.txt'.format(storage_prefix)
        if os.path.isfile(commands_filename):
            commands = []
            with open(commands_filename):
                for line in infile:
                    commands.append(line.strip().split(' '))
        return commands
    
    def run_all_commands(self, commands):
        processes = {}
        process_count = 1
        
        # initialize jobs
        for i in range(1, self.system.cores+1):
            processes[i] = subprocess.Popen(command.pop()+[str(i)]) # we must send unique running number
            
        # constantly check and make sure we are running the appropriate number of jobs
        while commands:
            for k, proc in processes.items():
                if proc.poll() != None:
                    processes[k] = subprocess.Popen(command+[str(i)])
            time.sleep(30)
            
        return