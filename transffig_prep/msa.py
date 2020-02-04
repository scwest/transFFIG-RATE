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
        self.commands_filename = ''
        
    def check_previous(self, commands, storage_prefix):
        self.commands_filename = '{}commands.txt'.format(storage_prefix)
        if os.path.isfile(self.commands_filename):
            commands = []
            with open(self.commands_filename):
                for line in infile:
                    commands.append(line.strip().split(' '))
        return commands
    
    def run_all_commands(self, commands):
        processes = {}
        process_count = 1
        
        # initialize jobs
        for i in range(1, self.system.cores+1):
            if commands:
                command = commands.pop()
                print(' '.join(command))
                processes[i] = subprocess.Popen(command+[str(i)]) # we must send unique running number
                raise Exception # testing one at a time
            
        # constantly check and make sure we are running the appropriate number of jobs
        while commands:
            for k, proc in processes.items():
                if proc.poll() != None:
                    command = commands.pop()
                    print(' '.join(command))
                    processes[k] = subprocess.Popen(commands+[str(i)])
            
            # keep record of unfinished commands
            self.update_command_file(commands)
            
            time.sleep(30)
            
            
        return
    
    def update_command_file(self, commands):
        with open(self.commands_filename, 'w') as outfile:
            for command in commands:
                outfile.write(' '.join(command) +'\n')
        return