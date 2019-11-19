# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 09:39:39 2019

@author: tombolelli
"""

from ase.io import read, Trajectory, write

#import traj
traj = Trajectory('tip3p_729mol_equil.traj')
#set watbox to last point of traj
watbox = traj[-1]
#get cell of watbox
Wcell = watbox.get_cell()

#import cluster coordinates
qm = read('coord')
#set cluster cell to watbox value
qm.set_cell(Wcell)
#center cluster in the cell
qm.center(vacuum=Wcell[0][0]/3)

#merge watbox and cluster
all = watbox.copy()
all.extend(qm)
all.write('temp.xyz')

#define index of atoms in mm part and index of atoms in qm part
mm_idx = list(range(0,(len(all) - len(qm))))
qm_idx = list(range((len(all) - len(qm)),len(all)))

#get distances between atoms in mm and atoms in qm
#get index of atoms within 2Ang from any qm atom
dist = []
conn_list = []

for i in mm_idx:
    for j in qm_idx:
        dist.append(all.get_distance(i, j))
        
for i, distance in enumerate(dist):
    i = int(i/len(qm_idx))
    if distance <=2:
        conn_list.append(i)

conn_set = set(conn_list)
conn_list = list(conn_set)

#get index of oxygens within 2Ang from any qm atom
oindex = []

for i in conn_list:
    if all[i].symbol == 'O':
        oindex.append(i)
        
#get index of water molecules within 2Ang from any qm qtom
watslist = []
for oxy in oindex:
    watslist.append(oxy)
    watslist.append(oxy +1)
    watslist.append(oxy +2)

watsset = set(watslist)
watslist = list(watsset)

#delete water molecules within 2Ang from any qm atom from all
del all.constraints
del all[watslist]

#write all to a file
all.write('wat-zn.xyz')
#reset constraints on waters!




