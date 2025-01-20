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
import sys
import shutil
def main():
    foldername = sys.argv[1]
    newfolder = sys.argv[2]
    previousfile = sys.argv[3]
    if previousfile!= 'none':
        oldf = open(previousfile, 'r')
        oldfilecontent = oldf.readlines()
        oldf.close()
    if os.path.isdir(newfolder) ==False:
        os.mkdir(newfolder)
    i = 1
    for ff in os.listdir(foldername):
        fullfolder = foldername + '/' + ff
        f = open(fullfolder +'/path.txt','r') #open the path, find successes
        content = f.readlines()
        f.close()
        for line in content:
            if 'success' in line:
                identification = line.split(' ')[2].strip(',')
                shutil.copy(fullfolder + '/' +identification + '.cpt', newfolder + '/' + str(i) + '.cpt' )
                if previousfile != 'none':
                    oldname = identification.split('_')[5]
                    for LLine in oldfilecontent:
                        if LLine.split(' , ')[-1].strip() == oldname.strip():
                            previousEntry = LLine[:-2] # this is to avoid a double \n situation, also delete the space
                f = open(newfolder +'/gronames.txt', 'a' )
                if previousfile != 'none':
                    f.write(previousEntry + ' , ' + identification + ' , ' + str(i) +' \n')
                else:
                    f.write(identification + ' , ' + str(i) +' \n')
                f.close()
                i +=1
    os.rename('controller_outputs.txt', newfolder + '/controller_outputs.txt')

    path = open(newfolder + 'path.txt', 'w')
    for ff in os.listdir(foldername):
        if os.path.isdir(foldername + '/' + ff):
            fullpath = foldername + '/' + ff
            for fil in os.listdir(fullpath):
                if fil[0] == 'C':
                    os.rename(fullpath + '/' + fil, newfolder +'/'+ fil)
                elif fil =='path.txt':
                    p1 = open(fullpath +'/' + fil, 'r')
                    p1c = p1.read()
                    p1.close()
                    path.write(p1c)
                elif fil == 'iterationtracker.txt':
                    os.rename(fullpath + '/' + fil, newfolder + '/folder_' + ff+ '_' +fil)
                elif fil =='outputfiles.txt':
                    os.rename(fullpath + '/' + fil, newfolder + '/folder_' + ff+ '_' +fil)
                elif fil == 'seedvalues.txt':
                    os.rename(fullpath + '/' + fil, newfolder + '/folder_' + ff+ '_' +fil)
if __name__ == "__main__":
    main()

