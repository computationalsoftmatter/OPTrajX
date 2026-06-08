#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH -J REX1basin_1
#SBATCH -t 280:00:00
#SBATCH --mail-user=manuel.lopez@yale.edu
#SBATCH --mail-type=end
#SBATCH -p pi_balou
# instead of pi_balou, we can also use the scavenge queue as needed.
# module restore
#module load GROMACS/2016.5-foss-2018b

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


import numpy as np
import os
import random
import sys
import subprocess
def replace_line(foldername, file_name, line_num, nofiles):
    if file_name not in os.listdir(foldername):
        out = open(foldername + '/' +file_name, 'w')
        for i in range(nofiles):
            if i == line_num:
                out.write('1 \n')
            else:
                out.write('0 \n')
        output = 1
    else:
        lines = open(foldername + '/' +file_name, 'r')
        content = lines.readlines()
        output = int(content[line_num]) +1
        content[line_num] = str(int(content[line_num]) +1) +'\n'
        out = open(foldername + '/' + file_name, 'w')
        out.writelines(content)
    out.close()
    return output
def PLMD(foldername, outputname): 
    line = -1
    f = open("plumed.dat", "r")
    content = f.readlines()
    s = content[line]
    out = ''  
    for i in range(len(s.split(' '))):
        if i!= 2:
            out += s.split(' ')[i] + ' '
        if i== 2:
            out += s.split(' ')[i].split('=')[0] + '=' + foldername+ '/C' + outputname + ' '
    content[line] = out
    outputfile = open(foldername+'/' + outputname+'.dat', 'w')
    outputfile.writelines(content)
    outputfile.close()
    f.close() 
def copier(original,later):
    if os.path.isfile(later)==False:
        lines = open( original, 'r')
        content = lines.readlines()
        lines.close()
        out = open(later, 'w')
        out.writelines(content)
        out.close()
def main():
    if len(sys.argv) ==9 or len(sys.argv) ==10: 
        for iteration in range(int(sys.argv[3])):
            foldername1 = sys.argv[1] 
            foldername2 = sys.argv[2]
            leftbound = float(sys.argv[4])
            rightbound = float(sys.argv[5])
            filepath = sys.argv[6]
            foldername = foldername1 + '/' +  foldername2
            if os.path.isdir(foldername1) == False:
                os.mkdir(foldername1)
            if iteration ==0:
                if os.path.isdir(foldername) == False: #this shouldn't ever fire but just in case
                    os.mkdir(foldername)
            f = open(foldername + "/outputfiles.txt", "a") #have moved this to smaller folder

            lis2 = os.listdir(filepath)
            id =  '.' + sys.argv[7] #".gro"    id is the keyword used for grouping runs, in this draft I just made .gro
            lis = [i for i in lis2 if id in i and 'bck' not in i and 'max' not in i]
            random.seed('folder'+ sys.argv[1]+ 'slurm' +str(sys.argv[2]) +'iteration' + str(iteration) )
            selected = random.choice(lis)
            if foldername+ "/grofiles.txt" not in os.listdir(): #make a list of .gro files to help with indexing in post processing
                grofiles = open(foldername+"/grofiles.txt", "w")
                for l in lis:
                    grofiles.write(str(l) + '\n')
                grofiles.close()

            number = replace_line(foldername,'iterationtracker.txt', lis.index(selected), len(lis))
            if len(selected.rsplit('.')) ==2:
                newname = "folder_" + sys.argv[1] + "_subfolder_"+sys.argv[2] + "_initialconfiguration_"+ selected.rsplit('.')[0] + '_iteration_' + str(number) #nohashes 
            else:
                newname = "folder_" + sys.argv[1] + "_subfolder_"+sys.argv[2] + "_initialconfiguration_"+ selected.rsplit('.')[0] +'.' + selected.rsplit('.')[1] + '_iteration_' + str(number)
            fullpath = os.getcwd() + '/'
            PLMD(fullpath+foldername, newname ) #edit plumed
            
            f.write(selected + '\n')
            f.close()

            lines = open( 'ffs.mdp', 'r') #mdp editer
            content = lines.readlines()
            lines.close()
            
            a = random.random()
            a = int(a*10**(len(str(a))-2))
            seedvalues = open(foldername + '/seedvalues.txt', 'a') #all txt files should go to big folder
            seedvalues.write(str(a) + '\n')
            seedvalues.close()
            content[-7] = content[-7].split('= ')[0] + '= ' + str(a) +'\n'
            out = open(foldername+'/'+newname+'.mdp', 'w')
            out.writelines(content)
            out.close() #end of mdp editer
            
            
            f2 = open(foldername +'/path.txt', "a")
            f2.write(selected[:-4] + ', '+ str(number)+', ' + newname + ' \n' ) #path file maker
            f2.close()
            id = "bck"    #after the run is complete delete all slurm and bck files
            lis2 = os.listdir()
            lis = [i for i in lis2 if id in i]
            for l in lis:
                os.system("rm "+ l)
            os.system('cd ' + fullpath+foldername+'/') #commands
            if sys.argv[7]== "gro":
                variabletopology = sys.argv[8]
                if variabletopology == '0': #this might not work, unclear if will be sent as string. I think it will
                    os.system("gmx_mpi grompp -f " +fullpath +foldername+'/'+ newname+'.mdp ' + "-c " +filepath +'/' + selected+ " -r " +filepath +'/'+ selected + " -p " + fullpath+ 'topol.top -n ' + fullpath +"index.ndx -o "   +fullpath +foldername+'/'+  newname+'.tpr ' + "-maxwarn 2")
                else:
                    os.system("gmx_mpi grompp -f " +fullpath +foldername+'/'+ newname+'.mdp ' + "-c " +filepath +'/' + selected+ " -r " +filepath +'/'+ selected + " -p " + variabletopology+'/'+ selected[:-4] + '.top -n ' + variabletopology +'/'+ selected[:-4] +".ndx -o "   +fullpath +foldername+'/'+  newname+'.tpr ' + "-maxwarn 2")
                os.system("srun gmx_mpi mdrun -v -deffnm " + fullpath +foldername+'/'+newname + " -plumed " +fullpath +foldername+'/'+ newname+'.dat -append -cpt 5' )
            if sys.argv[7] == "cpt":
                grofile = sys.argv[8]
                variabletopology = sys.argv[9]
                if variabletopology == '0':
                    os.system("gmx_mpi grompp -f "+fullpath + foldername + '/' + newname + '.mdp -c ' + grofile + ' -t ' + filepath + '/' + selected + " -r " +grofile+ " -p " + fullpath+ 'topol.top -n ' + fullpath +"index.ndx -o " + fullpath +foldername+'/'+  newname+'.tpr ' + "-maxwarn 2")
                else:
                    p = subprocess.run(["cat", filepath + "/gronames.txt"], capture_output=True, text=True)
                    a =p.stdout.strip().splitlines()
                    for i in range(len(a)):
                        if i != len(a)-1:
                            if a[i].split(' , ')[-1][:-1] ==selected[:-4]: #delete the .cpt also delete the space from a
                                firstname = a[i].split(' , ')[0]
                                topolname = firstname.split('initialconfiguration_')[1].split('_iteration')[0] #isolate the first groname
                        else: 
                            if a[i].split(' , ')[-1] ==selected[:-4]: #delete the .cpt also delete the space from a
                                firstname = a[i].split(' , ')[0]
                                topolname = firstname.split('initialconfiguration_')[1].split('_iteration')[0] #isolate the first groname
                    os.system("gmx_mpi grompp -f "+fullpath + foldername + '/' + newname + '.mdp -c ' + grofile + ' -t ' + filepath + '/' + selected + " -r " +grofile+ " -p " + variabletopology +'/'+ topolname + '.top -n ' + variabletopology + '/' + topolname + ".ndx -o " + fullpath +foldername+'/'+  newname+'.tpr ' + "-maxwarn 2")
                os.system("srun gmx_mpi mdrun -v -deffnm " + fullpath +foldername+'/'+newname + " -plumed " +fullpath +foldername+'/'+ newname+'.dat' )
                    


            if os.path.isfile(foldername + '/C'+ newname):#make sure this matches PLMD editer
                f = open(foldername + '/C'+ newname,'r' ) #let's see if this works
                content = f.readlines()
                f.close()
                f3 = open(foldername+'/path.txt', 'r') #found in shallow folder
                ct = f3.readlines()
                f3.close()
                c = content[-1].split(' ')[2]
                if float(c) >= rightbound: #success/ failure
                    ct[iteration] = ct[iteration].split(newname)[0] +newname+ ', success ' + ' \n'
                elif float(c)<= leftbound:
                    ct[iteration] = ct[iteration].split(newname)[0] +newname+ ', failure ' + ' \n'
                    #lis2 = os.listdir(foldername)
                    #lis = [i for i in lis2 if newname in i] #delete failed files
                    #for l in lis:
                        #os.system("rm "+foldername + '/' +  l) #although this code was helpful in clearing up data I think it will be more effective to delete in post to retain info on average length of failures
                else:
                    ct[iteration] = ct[iteration].split(newname)[0] +newname+ ', incomplete ' + ' \n'
                f4 = open(foldername +'/path.txt', 'w')
                f4.writelines(ct)
                f4.close()
            else:
                f3 = open(foldername+'/path.txt', 'r')
                ct = f3.readlines()
                f3.close()
                ct[iteration] = ct[iteration].split(newname)[0] +newname+ ', unclear' + ' \n'
                f4 = open(foldername+'/path.txt', 'w')
                f4.writelines(ct)
                f4.close()
if __name__ == "__main__":
    main()
