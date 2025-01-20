#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH -J REX1basin_1
#SBATCH -t 280:00:00
#SBATCH --mail-user=manuel.lopez@yale.edu
#SBATCH --mail-type=end
#SBATCH -p pi_balou

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

import os
import sys
import random
import subprocess
def runner2editor(foldername, outputname):
    line = 1
    f = open("runner2.py", "r")
    content = f.readlines()
    content[line] = content[line].split('=')[0] + '=\'' + outputname +'\' \n'

    outputfile = open(foldername+ '/' + outputname + '.py', 'w')
    outputfile.writelines(content)
    outputfile.close()
    f.close()
def PLMD(foldername, outputname): #finish cmd editor and make it so it's edited only in the folder
    line = -1
    f = open("plumed.dat", "r")
    content = f.readlines()
    s = content[line]
    out = ''  
    for i in range(len(s.split(' '))):
        if i!= 2:
            out += s.split(' ')[i] + ' '
        if i== 2:
           out += s.split(' ')[i].split('=')[0] + '=' + foldername+ '/RC' + outputname + ' ' #edit so you don't replace colvars
    content[line] = out
    outputfile = open(foldername+'/' + outputname+'.dat', 'w')
    outputfile.writelines(content)
    outputfile.close()
    f.close() #fix mdp, right now it is ffs.mdp 
def copier(original,later):
    if os.path.isfile(later)==False:
        lines = open( original, 'r')
        content = lines.readlines()
        lines.close()
        out = open(later, 'w')
        out.writelines(content)
        out.close()
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

def PLMDrunner(foldername, outputname): 
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

def main():
    if len(sys.argv) ==9 or len(sys.argv) ==10: 
        if os.path.isdir(sys.argv[1] +'/'+ sys.argv[2])==False:
            os.system('python3 runner.py ' + sys.argv[1] + ' ' + sys.argv[2] + ' ' + sys.argv[3] + ' ' + sys.argv[4] + ' ' + sys.argv[5])
        fullpath = os.getcwd() +'/'
        foldername1 = sys.argv[1]
        foldername2 = sys.argv[2]
        leftbound = float(sys.argv[4])
        rightbound = float(sys.argv[5])
        filepath = sys.argv[6]
        print(sys.argv[3])
        foldername = foldername1 + '/'+ foldername2


        f= open(foldername + '/path.txt', 'r') #in case parent didn't finish (write if it was a success/failure)
        content = f.readlines()
        f.close()
        it = -1
        for line in content:
            it = it +1
            le = line.split(' ')[2]
            if line.split(le)[1] ==' \n':
                if os.path.isfile(foldername + '/C'+ le): #if colvar file exists
                    f = open(foldername +'/C' + le,'r')
                    content = f.readlines()
                    f.close()
                    f3 = open(foldername+'/path.txt', 'r') #found in shallow folder
                    ct = f3.readlines()
                    f3.close()
                    c = content[-1].split(' ')[2]
                    if float(c) >= rightbound: #success/ failure
                        ct[it] = ct[it].split(le)[0] +le+ ', success ' + ' \n'
                    elif float(c)<= leftbound:
                        ct[it] = ct[it].split(le)[0] +le+ ', failure ' + ' \n'
                        #lis2 = os.listdir(foldername)
                        #lis = [i for i in lis2 if LastEntry in i] #delete failed files
                        #for l in lis:
                        #    os.system("rm "+foldername + '/' +  l)
                    else:
                        ct[it] = ct[it].split(le)[0] +le+ ', incomplete ' + ' \n'
                    f4 = open(foldername +'/path.txt', 'w')
                    f4.writelines(ct)
                    f4.close()
                else:
                    if os.path.isfile(foldername + '/' + le +'.gro'):
                        f3 = open(foldername+'/path.txt', 'r')
                        ct = f3.readlines()
                        f3.close()
                        ct[it] = ct[it].split(le)[0] +le+ ', unclear ' + ' \n'
                        f4 = open(foldername+'/path.txt', 'w')
                        f4.writelines(ct)
                        f4.close()
                    else:
                        f3 = open(foldername+'/path.txt', 'r')
                        ct = f3.readlines()
                        f3.close()
                        ct[it] = ct[it].split(le)[0] +le+ ', incomplete ' + ' \n'
                        f4 = open(foldername+'/path.txt', 'w')
                        f4.writelines(ct)
                        f4.close()


        for fil in os.listdir(foldername): #in case restarter didn't finish
            if fil[0]=='R':
                if os.path.isfile(foldername +'/' +fil[1:]): #check if C... exists
                    LastEntry = fil.split('RC')[1]
                    file1 = open(foldername + '/C'+ LastEntry,'r')
                    file2 = open(foldername + '/RC' + LastEntry,'r')
                    content1 = file1.readlines()
                    content2 = file2.readlines()
                    for i in range(1,len(content2)):
                        content2[i] = content2[i].split(' ')[0] + ' ' + str( float( content2[i].split(' ')[1]) + float(content1[-1].split(' ')[1] ) ) + ' ' + content2[i].split(' ')[2]
                    file1.close()
                    file2.close()
                    content3 = content1 + content2[2:]
                    os.system('rm ' + foldername + '/C' + LastEntry)
                    os.system('rm ' + foldername + '/RC' + LastEntry)
                    newfile = open(foldername + '/C'+ LastEntry,'w') #overide the colvar file
                    newfile.writelines(content3)
                    newfile.close()

                    f = open(foldername + '/C'+ LastEntry,'r' ) #let's see if this works
                    content = f.readlines()
                    f.close()
                    f3 = open(foldername+'/path.txt', 'r') #found in shallow folder
                    ct = f3.readlines()
                    f3.close()
                    c = content[-1].split(' ')[2] #assume if restarter didn't work then it is the last option in path that was restarted
                    if float(c) >= rightbound: #success/ failure
                        ct[-1] = ct[-1].split(LastEntry)[0] +LastEntry+ ', success ' + ' \n'
                    elif float(c)<= leftbound:
                        ct[-1] = ct[-1].split(LastEntry)[0] +LastEntry+ ', failure ' + ' \n'
                        #lis2 = os.listdir(foldername)
                        #lis = [i for i in lis2 if LastEntry in i] #delete failed files
                        #for l in lis:
                        #    os.system("rm "+foldername + '/' +  l)
                    else:
                        ct[-1] = ct[-1].split(LastEntry)[0] +LastEntry+ ', incomplete ' + ' \n'
                    f4 = open(foldername +'/path.txt', 'w')
                    f4.writelines(ct)
                    f4.close()
                else: #if the C file doesn't exist rename the r file and then follow the procedure for checking an unfinished parent file
                    os.rename(foldername + '/' + fil, foldername + '/' + fil[1:])
                    if os.path.isfile(foldername + '/' + LastEntry +'.gro'):
                        f3 = open(foldername+'/path.txt', 'r')
                        ct = f3.readlines()
                        f3.close()
                        ct[-1] = ct[-1].split(LastEntry)[0] +LastEntry+ ', unclear ' + ' \n'
                        f4 = open(foldername+'/path.txt', 'w')
                        f4.writelines(ct)
                        f4.close()
                    else:
                        f3 = open(foldername+'/path.txt', 'r')
                        ct = f3.readlines()
                        f3.close()
                        ct[-1] = ct[-1].split(LastEntry)[0] +LastEntry+ ', incomplete ' + ' \n'
                        f4 = open(foldername+'/path.txt', 'w')
                        f4.writelines(ct)
                        f4.close()


        f = open(foldername + '/path.txt')
        content = f.readlines()
        f.close()
        pathlength = len(content) 
        print(pathlength)
        print(int(sys.argv[3]))

            
        if pathlength< int(sys.argv[3]):
            for iteration in range(pathlength,int(sys.argv[3])): #this is the entire runner.py script
                foldername = foldername1 + '/' +  foldername2
                if os.path.isdir(foldername1) == False: # this section is in case the original runner.py didn't develop the demanded number of trajectories, will run it again for that many trajectories
                    os.mkdir(foldername1)
                if iteration ==0:
                    if os.path.isdir(foldername) == False: #this shouldn't ever fire but just in case
                        os.mkdir(foldername)
                f = open(foldername + "/outputfiles.txt", "a") #have moved this to smaller folder

                lis2 = os.listdir(filepath)
                id =  '.' + sys.argv[7] #".gro"    #id is the keyword used for grouping runs, in this draft I just made .gro
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
                PLMDrunner(fullpath+foldername, newname ) #edit plumed
                
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
                if  len(selected.rsplit('.')) ==2:
                    f2.write(selected.rsplit('.')[0] + ', '+ str(number)+', ' + newname + ' \n' ) #path file maker
                else:
                    f2.write(selected.rsplit('.')[0]  +'.' + selected.rsplit('.')[1] + ', '+ str(number)+', ' + newname + ' \n' ) #path file maker
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
                
            
        f = open(foldername + '/path.txt')
        content = f.readlines()
        f.close()
        pathlength = len(content) # this section is in case the original runner.py didn't develop the demanded number of trajectories, will run it again for that many trajectories
        print(pathlength)
        print(int(sys.argv[3]))
        for iteration in range(int(pathlength)):
            print(iteration)
            f3 = open(foldername +'/path.txt', 'r')  #need to get the LastEntry and selected
            fullcontent = f3.readlines()
            f3.close()
            LastEntry = fullcontent[iteration].split(' ')[2].strip(',')
            selected = fullcontent[iteration].split(' ')[0].strip(',') #do not update the restarter file name
            follower = fullcontent[iteration].split(LastEntry)[1].split(' ')[0]+ ' ' + fullcontent[iteration].split(LastEntry)[1].split(' ')[1]
            # check path for LastEntry, if words following it are success or failure, we will not run this code.
            if follower != ', success' and follower != ', failure':
                if os.path.isfile(foldername +  '/' +LastEntry + '.cpt'): #if cpt file exists...
                    f2 = open(foldername + '/' + LastEntry + '.mdp' , 'r') #open mdp file
                    content = f2.readlines()
                    content[-9] = content[-9].split('=')[0] + '= ' + 'no \n' #edit these
                    if content[-8][0]!= ';':
                        content[-8] = ';' + content[-8]
                    if content[-7][0]!= ';':
                        content[-7] = ';' + content[-7]
                    f4 = open(foldername + '/' + LastEntry + '.mdp' , 'w') #we are replacing the mdp file
                    f4.writelines(content)
                    f2.close()
                    f4.close() #this edits the mdp file and puts it in the correct folder
                    #deleted the part that edits the cmd file
                    PLMD(fullpath + foldername, LastEntry)  #this replaces the plmd file
                    #os.system('rm ' + foldername +'/'+  LastEntry + '.edr') #this needs to be combined
                    #if os.path.isfile(foldername + '/' + LastEntry +'.gro'):
                    #    os.system('rm ' + foldername +'/'+  LastEntry + '.gro') #this should be deleted
                    os.system('rm ' + foldername +'/'+  LastEntry + '.log') #log should be deleted
                    os.system('rm ' + foldername +'/'+  LastEntry + '.tpr') #tpr should be deleted
                    #os.system('rm ' + foldername +'/'+  LastEntry + '.trr') #trr will be combined
                    #os.system('rm ' + foldername +'/'+ LastEntry + '.xtc') #xtc will not occur in latest version of code
                    id = "bck"    #after the run is complete delete all slurm and bck files
                    lis2 = os.listdir()
                    lis = [i for i in lis2 if id in i]
                    for l in lis:
                        os.system("rm "+ l)
                    if sys.argv[7]== "gro":
                        variabletopology = sys.argv[8]
                        if variabletopology == '0':
                            entry = ('gmx_mpi grompp -f ' + fullpath + foldername + '/' + LastEntry + '.mdp -c ' + filepath+'/' + selected + '.gro -t ' + fullpath + foldername + '/' + LastEntry + '.cpt ' + '-r ' + filepath +'/' + selected + '.gro -p ' + fullpath + 'topol.top -n ' + fullpath + 'index.ndx -o ' + fullpath + foldername + '/' + LastEntry + '.tpr -maxwarn 2')
                        else:
                            entry = ('gmx_mpi grompp -f ' + fullpath + foldername + '/' + LastEntry + '.mdp -c ' + filepath+'/' + selected + '.gro -t ' + fullpath + foldername + '/' + LastEntry + '.cpt ' + '-r ' + filepath +'/' + selected + '.gro -p ' + variabletopology+'/'+ selected[:-4] + '.top -n ' + variabletopology +'/'+ selected[:-4] +".ndx -o " + fullpath + foldername + '/' + LastEntry + '.tpr -maxwarn 2')
                        #print(entry)
                        os.system(entry)
                        #print("srun gmx_mpi mdrun -v -deffnm " + fullpath +foldername+'/'+LastEntry + " -plumed " +fullpath +foldername+'/'+ LastEntry+'.dat')
                        os.system("srun gmx_mpi mdrun -v -deffnm " + fullpath +foldername+'/'+LastEntry + " -plumed " +fullpath +foldername+'/'+ LastEntry+'.dat' )
                    if sys.argv[7] == "cpt":
                        grofile = sys.argv[8]
                        variabletopology = sys.argv[9]
                        if variabletopology == '0':
                            os.system("gmx_mpi grompp -f "+fullpath + foldername + '/' + LastEntry + '.mdp -c ' + grofile + ' -t ' + filepath + '/' + selected + " -r " +grofile+ " -p " + fullpath+ 'topol.top -n ' + fullpath +"index.ndx -o " + fullpath +foldername+'/'+  LastEntry+'.tpr ' + "-maxwarn 2")
                        else:
                            p = subprocess.run(["cat", filepath + "/gronames.txt"], capture_output=True, text=True)
                            a =p.stdout.strip().splitlines()
                            for i in range(len(a)):
                                if i != len(a) -1:
                                    if a[i].split(' , ')[-1][:-1] ==selected[:-4]: #delete the .cpt also delete the space from a
                                        firstname = a[i].split(' , ')[0]
                                        topolname = firstname.split('initialconfiguration_')[1].split('_iteration')[0] #isolate the first groname
                                else:
                                    if a[i].split(' , ')[-1] ==selected[:-4]: #delete the .cpt also delete the space from a
                                        firstname = a[i].split(' , ')[0]
                                        topolname = firstname.split('initialconfiguration_')[1].split('_iteration')[0] #isolate the first groname
                            os.system("gmx_mpi grompp -f "+fullpath + foldername + '/' + LastEntry + '.mdp -c ' + grofile + ' -t ' + filepath + '/' + selected + " -r " +grofile+ " -p " +variabletopology +'/'+ topolname+ '.top -n ' + variabletopology +'/'+ topolname +".ndx -o " + fullpath +foldername+'/'+  LastEntry+'.tpr ' + "-maxwarn 2")
                        os.system("srun gmx_mpi mdrun -v -deffnm " + fullpath +foldername+'/'+LastEntry + " -plumed " +fullpath +foldername+'/'+ LastEntry+'.dat' )
                    
                else: #if cpt doesn't exist then we should just trigger the cmd file with new outputs and edit the path accordingly
                    if os.path.isfile(foldername +  '/C' +LastEntry ):
                        PLMD(fullpath + foldername, LastEntry)  #unsure if this is even necessary, would there be colvar outputs if there aren't cpt outputs???
                    if sys.argv[7]== "gro":
                        variabletopology = sys.argv[8]
                        if variabletopology == '0':
                            os.system("gmx_mpi grompp -f " +fullpath +foldername+'/'+ LastEntry+'.mdp ' + "-c " +filepath +'/' + selected+ '.gro' + " -r " +filepath +'/'+ selected +'.gro'+  " -p " + fullpath+ 'topol.top -n ' + fullpath +"index.ndx -o "   +fullpath +foldername+'/'+ LastEntry+'.tpr ' + "-maxwarn 2")
                        else:
                            os.system("gmx_mpi grompp -f " +fullpath +foldername+'/'+ LastEntry+'.mdp ' + "-c " +filepath +'/' + selected+ '.gro' + " -r " +filepath +'/'+ selected +'.gro'+  " -p " + variabletopology+'/'+ selected[:-4] + '.top -n ' + variabletopology +'/'+ selected[:-4] +".ndx -o "  +fullpath +foldername+'/'+ LastEntry+'.tpr ' + "-maxwarn 2")
                        os.system("srun gmx_mpi mdrun -v -deffnm " + fullpath +foldername+'/'+LastEntry + " -plumed " +fullpath +foldername+'/'+ LastEntry+'.dat' )
                    if sys.argv[7] == "cpt":
                        grofile = sys.argv[8]
                        variabletopology = sys.argv[9]
                        if variabletopology == '0':
                            os.system("gmx_mpi grompp -f "+fullpath + foldername + '/' + LastEntry + '.mdp -c ' + grofile + ' -t ' + filepath + '/' + selected + " -r " +grofile+ " -p " + fullpath+ 'topol.top -n ' + fullpath +"index.ndx -o " + fullpath +foldername+'/'+  LastEntry+'.tpr ' + "-maxwarn 2")
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
                            os.system("gmx_mpi grompp -f "+fullpath + foldername + '/' + LastEntry + '.mdp -c ' + grofile + ' -t ' + filepath + '/' + selected + " -r " +grofile+ " -p " +variabletopology +'/'+ topolname+ '.top -n ' + variabletopology +'/'+ topolname +".ndx -o "+ fullpath +foldername+'/'+  LastEntry+'.tpr ' + "-maxwarn 2")
                        os.system("srun gmx_mpi mdrun -v -deffnm " + fullpath +foldername+'/'+LastEntry + " -plumed " +fullpath +foldername+'/'+ LastEntry+'.dat' )

                #should add in something that merges he colvar files
                if os.path.isfile(foldername +  '/RC' +LastEntry ) and os.path.isfile(foldername +  '/C' +LastEntry ) : #if have multiple colvar files combine them
                    file1 = open(foldername + '/C'+ LastEntry,'r')
                    file2 = open(foldername + '/RC' + LastEntry,'r')
                    content1 = file1.readlines()
                    content2 = file2.readlines()
                    for i in range(1,len(content2)):
                        content2[i] = content2[i].split(' ')[0] + ' ' + str( float( content2[i].split(' ')[1]) + float(content1[-1].split(' ')[1] ) ) + ' ' + content2[i].split(' ')[2] 
                    file1.close()
                    file2.close()
                    content3 = content1 + content2[2:]
                    os.system('rm ' + foldername + '/C' + LastEntry)
                    os.system('rm ' + foldername + '/RC' + LastEntry)
                    newfile = open(foldername + '/C'+ LastEntry,'w') #overide the colvar file
                    newfile.writelines(content3)
                    newfile.close()
                if os.path.isfile(foldername + '/C'+ LastEntry):#make sure this matches PLMD editer
                    f = open(foldername + '/C'+ LastEntry,'r' ) #let's see if this works
                    content = f.readlines()
                    f.close()
                    f3 = open(foldername+'/path.txt', 'r') #found in shallow folder
                    ct = f3.readlines()
                    f3.close()
                    c = content[-1].split(' ')[2]
                    if float(c) >= rightbound: #success/ failure
                        ct[iteration] = ct[iteration].split(LastEntry)[0] +LastEntry+ ', success ' + ' \n'
                    elif float(c)<= leftbound:
                        ct[iteration] = ct[iteration].split(LastEntry)[0] +LastEntry+ ', failure ' + ' \n'
                    else:
                        ct[iteration] = ct[iteration].split(LastEntry)[0] +LastEntry+ ', incomplete ' + ' \n'
                    f4 = open(foldername +'/path.txt', 'w')
                    f4.writelines(ct)
                    f4.close()
                else:
                    f3 = open(foldername+'/path.txt', 'r')
                    ct = f3.readlines()
                    f3.close()
                    ct[iteration] = ct[iteration].split(LastEntry)[0] +LastEntry+ ', unclear ' + ' \n'
                    f4 = open(foldername+'/path.txt', 'w')
                    f4.writelines(ct)
                    f4.close()
                


if __name__ == "__main__":
    main()
