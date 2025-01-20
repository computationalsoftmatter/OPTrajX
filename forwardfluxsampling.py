"""
===========================================================
Forward Flux Sampling (FFS) Implementation
===========================================================

This Python code implements Forward Flux Sampling (FFS) for studying rare events
in stochastic systems. It includes modules for defining order parameters,
generating configurations, and computing transition probabilities.

License:
--------
This software is distributed under the GNU General Public License (GPL),
either version 3 of the License or (at your option) any later version.
You are free to use, modify, and distribute this software under the terms
of the GPL. For more details, see the LICENSE file included with this project
or visit https://www.gnu.org/licenses/gpl-3.0.html.

Contributors:
-------------
- Manuel Lopez
- Betul Uralcan
- Amir Haji-akbari

Contact:
--------
For questions, please contact:
- Betul Uralcan
  - Email: betul.uralcan@bogazici.edu.tr
  - Email: betul.uralcan@yale.edu

Usage:
------
This script is designed for educational and research purposes. Please cite
the appropriate references if you use this code in your research.

===========================================================
"""

import subprocess
import os
import time
import sys

def crossingcounter(foldername):
    noCrossings=0
    total = 0
    successlength = 0
    faillength = 0
    nofails = 0
    mx = 0
    
    paths = ''
    if os.path.isdir(foldername):
            for FF in os.listdir(foldername):
                if os.path.isdir(foldername +'/' + FF):
                    fullfolder = foldername +'/'+  FF
                    if os.path.isfile(fullfolder + '/path.txt'):
                        paths = paths + fullfolder+'/path.txt ' #first pipe the paths so that they aren't overwritten
    if paths != '':
        os.system('cat -t ' + paths + '>globalPath.txt')
        path = open('globalPath.txt', 'r')
        content = path.readlines()
        path.close()
        for line in content:
            total +=1
            if 'success' in line:
                noCrossings +=1
                identification = line.split(' ')[2].strip(',')
                subfolder = identification.split('_')[3] #with this piping mechanism only open colvar files that are complete
                if os.path.isfile(foldername + '/'+subfolder +'/C'+identification):
                    colvarf = open(foldername + '/'+subfolder +'/C'+identification,'r')
                    colvar = colvarf.readlines()
                    colvarf.close()
                    successlength += float(colvar[-1].split(' ')[1])
                    if float(colvar[-1].split(' ')[2])> mx:
                        mx = float(colvar[-1].split(' ')[2])
            if 'failure' in line:
                nofails +=1
                identification = line.split(' ')[2].strip(',')
                subfolder = identification.split('_')[3] #with this piping mechanism only open colvar files that are complete
                if os.path.isfile(foldername + '/'+subfolder +'/C'+identification):
                    colvarf = open(foldername +'/'+subfolder +'/C'+identification,'r')
                    colvar = colvarf.readlines()
                    colvarf.close()
                    faillength += float(colvar[-1].split(' ')[1])
        os.system('rm globalPath.txt') #this needs to be deleted so that it can be piped to later on
    return noCrossings,nofails,  total, successlength, faillength, mx

def main():
    entries = [0,0,0,0,0,0,0,0,0,0] #did it this way so I wouldn't have to import numpy (not a default package)
    filetypes = ['batch size', 'number of restarters', 'output folder name', 'trajectories per folder', 'desired number of crossings', 'basin order parameter', 'target order parameter', 'path to input files', 'run time', 'workflow (cpt or gro)' ]
    flags = ['-b', '-r', '-n', '-s', '-c', '-i', '-g', '-p', '-t', '-gro or -cpt']
    tt = '24:00:00'
    variabletopology = 0
    for i in range(len(sys.argv)): #flag system
        if sys.argv[i] =='-b': # b for batch size
            entries[0] = 1
            BatchSize =int( sys.argv[i+1])
        if sys.argv[i] =='-r': #r for restarts
            entries[1] = 1
            InitialSize =int( sys.argv[i+1])
        if sys.argv[i] =='-n': #n for name
            entries[2] = 1
            foldername =sys.argv[i+1]
        if sys.argv[i] =='-s': #s for size
            entries[3] = 1
            runsperfolder = int( sys.argv[i+1])
        if sys.argv[i] =='-c': #c for crossings
            entries[4] = 1
            desiredCrossings = int( sys.argv[i+1])
        if sys.argv[i] =='-i': #i for initial (basin)
            entries[5] = 1
            leftbound = sys.argv[i+1]
        if sys.argv[i] =='-g': #g for goal lambda
            entries[6] = 1
            rightbound = sys.argv[i+1]
        if sys.argv[i] =='-p': #p for path
            entries[7] = 1
            filespath = sys.argv[i+1]
        if sys.argv[i] == '-t': #t for time
            entries[8] = 1
            tt = sys.argv[i+1]
        if sys.argv[i] =='-gro': #for wants gro file workflow
            entries[9] = 1
            filetype = 'gro'
        if sys.argv[i] =='-cpt':
            entries[9] = 1
            filetype = 'cpt'
            grofile = sys.argv[i+1]
            if grofile.split('.')[-1] != 'gro' or os.path.isfile(grofile) ==False:
                print("Error, no valid gro file given for checkpoint workflow")
                sys.exit(1)
        if sys.argv[i] == '-x':
            variabletopology = sys.argv[i+1]

    listofmissingfiles = ''
    missingflags = ''
    no = 0
    for i in range(len(entries)):
        if entries[i] ==0 and i!=4 and i!=8:
            if no == 0:
                listofmissingfiles += filetypes[i]
                missingflags += flags[i]
            else:
                listofmissingfiles += ' and '+ filetypes[i]
                missingflags += ' and '+ flags[i] 
            no+=1
        
    if listofmissingfiles != '':
        print("Error, require input for " + listofmissingfiles)
        if no==1:
            print("Try inserting the following flag " + missingflags)
        else:
            print("Try inserting the following flags " + missingflags)
        sys.exit(1)

    if entries[4]==0:
        print('Warning proceeding without interupting once a certain number of crossings has been reached')
        desiredCrossings = BatchSize*runsperfolder
    if entries[8] ==0:
        print('Warning proceeding with default run time of 24 hours')
    if variabletopology ==0:
        print("Warning proceeding assuming constant topology file (topol.top) and index file (index.ndx), if this is not desired restart simulation using flag -x followed by the folder holding these files")
    


    joblist = []



    
    #adjust .dat file to have correct bounds
    f1 = open('plumed.dat', 'r')
    ct = f1.readlines()
    f1.close()
    line = 6
    ct[line] = ct[line].split('=')[0] + '=' + leftbound + '\n'
    ct[line+1] =ct[line+1].split('=')[0] + '=' + rightbound+ '\n'
    ct[line+2] =ct[line+2].split('=')[0] + '=' + str(float(rightbound) + .1)+ '\n'
    f2 = open('plumed.dat', 'w')
    f2.writelines(ct)
    f2.close()
    
    # Submitting parent_bash.sh job
    if filetype == "gro":
        parent_process = subprocess.run(["sbatch", "-N", "1", "-n", "28", "--time", tt, "-p", "scavenge","-a", "1-" +str(BatchSize), "parent_bash.sh", foldername,  str(runsperfolder), leftbound,rightbound, filespath, filetype,variabletopology], capture_output=True, text=True)
        s_old = parent_process.stdout.strip().split()[-1]
    elif filetype =="cpt":
        parent_process = subprocess.run(["sbatch", "-N", "1", "-n", "28", "--time", tt, "-p", "scavenge","-a", "1-" +str(BatchSize), "cpt_parent_bash.sh", foldername,  str(runsperfolder), leftbound,rightbound, filespath, filetype, grofile,variabletopology], capture_output=True, text=True)
        s_old = parent_process.stdout.strip().split()[-1]
    joblist.append(s_old)

    # Submitting restarter_bash.sh job


    for i in range(InitialSize):
        if filetype == "gro":
            restarter_process = subprocess.run(["sbatch", "-N", "1", "-n", "28", "--time", tt, "-p", "scavenge","-a","1-" + str(BatchSize), "restarter_bash.sh",foldername, str(runsperfolder),leftbound,rightbound,filespath,filetype,variabletopology], capture_output=True, text=True)
            s_new = restarter_process.stdout.strip().split()[-1]
        elif filetype =="cpt":
            restarter_process = subprocess.run(["sbatch", "-N", "1", "-n", "28", "--time", tt, "-p", "scavenge","-a","1-" + str(BatchSize), "cpt_restarter_bash.sh",foldername, str(runsperfolder),leftbound,rightbound,filespath,filetype,grofile,variabletopology], capture_output=True, text=True)
            s_new = restarter_process.stdout.strip().split()[-1]


    # Loop to update job dependencies
        for I in range(1, BatchSize+1):
            subprocess.run(["scontrol", "update", f"JobId={s_new}_{I}", f"Dependency=afterany:{s_old}_{I}"])
        s_old = s_new
        joblist.append(s_old)


    noCrossings = 0
    f = open("controller_outputs.txt", 'w')
    f.write("total, successful, percent, fails, percent, max λ,  ave success length, ave failure length \n"  ) #add in average length
    f.close()
    print("waiting for simulations to begin running")
    qC = 1
    while noCrossings< desiredCrossings and qC != 0:
        qC = 0
        queue = subprocess.run(['squeue', '--me'],  capture_output=True, text=True)
        q = queue.stdout #check to see if this worked
        for entry in joblist:
            if entry in q:
                qC =1

        id = "bck"    #after the run is complete delete all slurm and bck files
        lis2 = os.listdir()
        lis = [i for i in lis2 if id in i]
        for l in lis:
            os.system("rm "+ l)
        noCrossings,nofails,  total, successlength, faillength, mx  = crossingcounter(foldername)
        if noCrossings > 0:
            print('crossings= ' + str(noCrossings))
            if nofails ==0:
                f = open("controller_outputs.txt", "a")
                f.write(str(total) +' , ' + str(noCrossings) +' , '+ str(round(noCrossings/total, 4))+ ' , ' + str(nofails) + ' , ' + str(round(nofails/total,2)) + ' , '+ str(round(mx,2)) +' , ' + str(round(successlength/noCrossings,2)) + ' , ' + '0' +  ' \n' )
                f.close()
            else: 
                f = open("controller_outputs.txt", "a")
                f.write(str(total) +' , ' + str(noCrossings) +' , '+ str(round(noCrossings/total, 2))+ ' , ' + str(nofails) + ' , ' + str(round(nofails/total,2)) + ' , ' +str(round(mx,2)) +' , ' + str(round(successlength/noCrossings,2)) + ' , ' + str(round(faillength/nofails,2)) +  ' \n' )
                f.close()
        time.sleep(5)
    for job in joblist:
        for I in range(1,BatchSize+1):
            os.system("scancel " + job + '_' + str(I))
    id = "slurm"
    lis2 = os.listdir()
    lis = [i for i in lis2 if id in i]
    for l in lis:
        os.system("rm "+ l)
    id = "bck"    #after the run is complete delete all slurm and bck files
    lis2 = os.listdir()
    lis = [i for i in lis2 if id in i]
    for l in lis:
        os.system("rm "+ l)
    os.system('sbatch emailbash.sh')

if __name__ == "__main__":
    main()
