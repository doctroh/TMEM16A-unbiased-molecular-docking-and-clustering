import pubchempy as pcp
import os
import shutil
import paramiko
import time
from scp import SCPClient
from Bio.PDB import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import seaborn as sns


# create linux input folder
def linuxFolder (cid) :
    global receptor

    os.system ('if not exist "linux_%s" mkdir "linux_%s"' % (cid, cid))
    shutil.copy ("%s.pdbqt" % (receptor), "linux_%s" % (cid))
    return

# download ligand file
def getCoord (cid) :
    pcp.download ("SDF", "linux_%s/%s.sdf" % (cid, cid), cid, "cid", overwrite = True)
    return

# generate ligand pdbqt file
def SDFtoPDBQT (cid) :
    os.system ("obabel linux_%s/%s.sdf -O linux_%s/%s.pdbqt --gen3d" % (cid, cid, cid, cid))
    return

# generate conf file
def confWriteGlobal (cid) :
    global cpu, receptor

    coord = [0, 0, 0] # center of mass
    size = [150, 150, 150]
    cpu = 4
    exhaustiveness = 8
    num_modes = 10
    energy_range = 3

    f = open ("linux_%s/%s_conf.txt" % (cid, cid), "w")
    f.write ("receptor = %s.pdbqt\n" % (receptor))
    f.write ("ligand = %s.pdbqt\n" % (cid))
    f.write ("center_x = %f\n" % (coord [0]))
    f.write ("center_y = %f\n" % (coord [1]))
    f.write ("center_z = %f\n" % (coord [2]))
    f.write ("size_x = %f\n" % (size [0]))
    f.write ("size_y = %f\n" % (size [1]))
    f.write ("size_z = %f\n" % (size [2]))
    f.write ("cpu = %d\n" % (cpu))
    f.write ("exhaustiveness = %d\n" % (exhaustiveness))
    f.write ("num_modes = %d\n" % (num_modes))
    #f.write ("energy_range = %d\n" % (energy_range))
    f.close ()
    return

# generate linux input file
def linuxInput (cid) :
    global cpu

    node = "1"
    run = 10
    localRun = 10

    for K in range (1, run + 1) :
        #f = open ("linux_%s/%s_linuxRun_%d.sh" % (cid, cid, K), "w")
        f = open ("C:/Users/jwroh/Desktop/VM folder/%s_linuxRun_%d.sh" % (cid, K), "w")
        f.write ("#!/bin/tcsh\n")
        f.write ("#PBS -N Vina_%d\n" % (K))
        f.write ("#PBS -l nodes=%s:ppn=%d\n" % (node, cpu))
        f.write ("#PBS -q workq\n")
        f.write ("cd $PBS_O_WORKDIR\n")
        f.write ("\n")
        f.write ("mkdir -p logResult\n")
        f.write ("mkdir -p outResult\n")
        f.write ("\n")
        f.write ("set cid = %d\n" % (cid))
        f.write ("set cnt = 0\n")
        f.write ("set cntmax = %d\n" % (localRun - 1))
        f.write ("\n")
        f.write ("while (${cnt} <= ${cntmax})\n")
        f.write ("    vina --config %s_conf.txt --log logResult/%s_log_%d_${cnt}.txt --out outResult/%s_out_%d_${cnt}.pdbqt\n" % (cid, cid, K, cid, K))
        f.write ("    @ cnt += 1\n")
        f.write ("end\n")
        f.close ()

    f = open ("linux_%s/pyRun_%s.py" % (cid, cid), "w")
    f.write ("import os\n")
    f.write ("import time\n")
    f.write ("\n")
    for K in range (1, run + 1) :
        f.write ('os.system ("qsub %s_linuxRun_%d.sh")\n' % (cid, K))
        f.write ('time.sleep (1.0)\n')
    f.close ()
    return

'''
# connect to SSH server -> abandoned due to instability and unavailability of dos2unix program
def SSH () :
    global cid, run

    cli = paramiko.SSHClient ()
    cli.set_missing_host_key_policy (paramiko.AutoAddPolicy)

    server = "203.252.74.112"
    user = "jwroh" 
    pwd = "tmem16"

    cli.connect (server, port = 22, username = user, password = pwd)
    channel = cli.invoke_shell ()
    scp = SCPClient(cli.get_transport())

    for file in os.listdir ("linux_%s" % (cid)):
        scp.put ("linux_%s/%s" % (cid, file), "temp/")

    channel.send ("cd temp/\n")
    for i in range (1) :
        channel.send ("qsub %s_linuxRun_%d.sh\n" % (cid, i + 1))

    print ()
    output = channel.recv (65535).decode ("utf-8")
    print (output)

    cli.close ()
    return
SSH ()
'''






# create analysis output folder
def analysisFolder (cid) :
    os.system ('if not exist "analysis_%s" mkdir "analysis_%s"' % (cid, cid))
    return

# measure center of mass
def cenMass (cid) :
    f = open ("analysis_%s/cenMass.txt" % (cid), "w")
    for file in os.listdir ("linux_%s/outResult" % (cid)) :
        parser = PDBParser ()
        structure = parser.get_structure ("LIG", "linux_%s/outResult/%s" % (cid, file))
        for model in structure:
            for i in model.center_of_mass () :
                f.write ("%.3f\t" % (i))
            f.write ("\n")
    f.close ()
    return

# calculate based on C2 symmetry
def cenMassSym (cid) :
    global receptor

    parser = PDBParser ()
    structure = parser.get_structure (receptor, "%s.pdb" % (receptor))

    temp = []
    cen = []
    for model in structure :
        for i in model.center_of_mass () :
            cen.append (i)
        for chain in model :
            for residue in chain :
                if str (residue) [26 : 29] == "666" :
                    for i in residue.center_of_mass () :
                        temp.append (i)

    [cenX, cenY, cenZ] = cen
    [cenAX, cenAY, cenAZ] = temp [ : 3]
    [cenBX, cenBY, cenBZ] = temp [3 : ]

    f = open ("analysis_%s/cenMass.txt" % (cid), "r")
    g = open ("analysis_%s/cenMass_sym.txt" % (cid), "w")

    for line in f :
        line = line.strip ().split ("\t")
        temp = []
        for i in line :
            temp.append (float (i))
        [X, Y, Z] = temp
        
        disToCenA = (X - cenAX) ** 2 + (Y - cenAY) ** 2
        disToCenB = (X - cenBX) ** 2 + (Y - cenBY) ** 2

        if min (disToCenA, disToCenB) == disToCenA :
            g.write ("%.3f\t%.3f\t%.3f\t\n" % (X, Y, Z))
        elif min (disToCenA, disToCenB) == disToCenB :
            X1 = (X - cenX) * np.cos (np.deg2rad (180)) - (Y - cenY) * np.sin (np.deg2rad (180)) + cenX
            Y1 = (X - cenX) * np.sin (np.deg2rad (180)) + (Y - cenY) * np.cos (np.deg2rad (180)) + cenY
            g.write ("%.3f\t%.3f\t%.3f\t\n" % (X1, Y1, Z))

    f.close ()
    g.close ()
    return

# KMEANS cluster analysis
def KMEANScluster (cid) :
    point = []
    f = open ("analysis_%s/cenMass_sym.txt" % (cid), "r")
    for line in f :
        line = line.strip ().split ("\t")
        x = float (line [0])
        y = float (line [1])
        z = float (line [2])
        point.append ([x, y, z])
    f.close ()

    clusterNum = range (1, 20)
    iner = []

    cluster = 7
    while True :
        clus = KMeans (n_clusters = cluster).fit (point)
        iner.append (clus.inertia_)

        clusterLabel = clus.labels_
        clusterCount = np.bincount (clusterLabel)
        clusterCountSort = clusterCount.sort ()

        mean = []
        for i in range (cluster) :
            tempX = []
            tempY = []
            tempZ = []
            for k in range (len (clusterLabel)) :
                if clusterLabel [k] == i :
                    tempX.append (point [k] [0])
                    tempY.append (point [k] [1])
                    tempZ.append (point [k] [2])
            mean.append ([len (tempX), [np.mean (tempX), np.mean (tempY), np.mean (tempZ)]])

        count1 = 0
        count2 = 0
        dim = 12.5
        within = [0] * cluster

        for i in range (cluster) :
            tempX = []
            tempY = []
            tempZ = []
            for k in range (len (clusterLabel)) :
                if clusterLabel [k] == i :
                    tempX.append (point [k] [0])
                    tempY.append (point [k] [1])
                    tempZ.append (point [k] [2])

            count = 0
            if mean [i] [0] >= len (point) / cluster :
                count1 += 1
                for k in range (len (tempX)) :
                    if mean [i] [1] [0] - dim <= tempX [k] <= mean [i] [1] [0] + dim :
                        if mean [i] [1] [1] - dim <= tempY [k] <= mean [i] [1] [1] + dim :
                            if mean [i] [1] [2] - dim <= tempZ [k] <= mean [i] [1] [2] + dim :
                                count += 1
                temp = float (count / mean [i] [0])
                if temp >= 0.9 :
                    count2 += 1
                    within [i] = temp

        if count1 == count2 :
            print ()
            print ("%s\t%d" % (cid, cluster))
            print (clusterCount)
            for i in range (cluster) :
                if within [i] != 0 :
                    print ("%d\t%.2f\t%.3f\t%.3f\t%.3f\t"% (mean [i] [0], within [i], mean [i] [1] [0], mean [i] [1] [1], mean [i] [1] [2]))
            print ()

            f = open ("analysis_%s/cenMass_sym_cluster.txt" % (cid), "w")
            for i in range (cluster) :
                if within [i] != 0 :
                    f.write ("%d\t%.2f\t%.3f\t%.3f\t%.3f\t\n"% (mean [i] [0], within [i], mean [i] [1] [0], mean [i] [1] [1], mean [i] [1] [2]))
            f.close ()
            break
        else :
            cluster += 1
    return

# generate conf file
def confWritePocket (cid) :
    global receptor

    size = [25, 25, 25]
    cpu = 8
    exhaustiveness = 8
    num_modes = 1
    energy_range = 3

    coord = []
    f = open ("analysis_%s/cenMass_sym_cluster.txt" % (cid), "r")
    for line in f :
        line = line.strip ().split ("\t")
        temp = []
        for i in line [2 : ] :
            temp.append (float (i))
        coord.append (temp)
    f.close ()

    shutil.copy ("%s.pdbqt" % (receptor), "analysis_%s" % (cid))
    shutil.copy ("linux_%s/%s.pdbqt" % (cid, cid), "analysis_%s" % (cid))

    for i in range (len (coord)) :
        f = open ("analysis_%s/pocketDocking_conf_%d.txt" % (cid, i), "w")
        f.write ("receptor = analysis_%s/%s.pdbqt\n" % (cid, receptor))
        f.write ("ligand = analysis_%s/%s.pdbqt\n" % (cid, cid))
        f.write ("center_x = %f\n" % (coord [i] [0]))
        f.write ("center_y = %f\n" % (coord [i] [1]))
        f.write ("center_z = %f\n" % (coord [i] [2]))
        f.write ("size_x = %f\n" % (size [0]))
        f.write ("size_y = %f\n" % (size [1]))
        f.write ("size_z = %f\n" % (size [2]))
        f.write ("cpu = %d\n" % (cpu))
        f.write ("exhaustiveness = %d\n" % (exhaustiveness))
        f.write ("num_modes = %d\n" % (num_modes))
        f.write ("energy_range = %d\n" % (energy_range))
        f.close ()
    return

# pocket specific docking
def pocketDocking (cid) :
    count = 0
    f = open ("analysis_%s/cenMass_sym_cluster.txt" % (cid), "r")
    for line in f :
        line = line.strip ().split ("\t")
        if len (line) != 0 :
            count += 1
    f.close ()

    for i in range (count) :
        print ("%s\t#%d" % (cid, i + 1))
        os.system ('"vina.exe" --config analysis_%s/pocketDocking_conf_%d.txt --log analysis_%s/pocketDocking_log_%d.txt --out analysis_%s/pocketDocking_out_%d.pdbqt' % (cid, i, cid, i, cid, i))
        print ()
    return

# docking result output
def resultOutput (cid) :
    result = []
    f = open ("analysis_%s/cenMass_sym_cluster.txt" % (cid), "r")
    for line in f :
        line = line.strip ().split ("\t")
        if len (line) != 0 :
            temp = []
            for i in line :
                temp.append (float (i))
            result.append (temp)
    f.close ()

    for i in range (len (result)) :
        f = open ("analysis_%s/pocketDocking_log_%d.txt" % (cid, i), "r")
        for line in f :
            line = line.strip ().split ()
            if len (line) != 0 and line [0] == "1" :
                result [i].append (float (line [1]))
        f.close ()

    result = sorted (result, key = lambda l:l [0], reverse = True)
    
    f = open ("analysis_%s/pocketDocking__final.txt" % (cid), "w")
    for i in range (len (result)) :
        f.write ("%.0f\t%.1f\t%.3f\t%.3f\t%.3f\t\n" % (result [i] [0], result [i] [5], result [i] [2], result [i] [3], result [i] [4]))
    f.close ()

    os.system ("obabel analysis_%s/pocketDocking_out_*.pdbqt -O analysis_%s/pocketDocking__final.pdb" % (cid, cid))
    return





receptor = "hT16A_m_dimer"
a = int (input ("Enter mode [1; input generation, 2; output analysis] : "))
CID = []
temp = input ("Enter PubChem CID : ").strip ().split ()
for i in temp :
    CID.append (int (i))

if a == 1 :
    for i in CID :
        linuxFolder (i)
        getCoord (i)
        SDFtoPDBQT (i)
        confWriteGlobal (i)
        linuxInput (i)

if a == 2 :
    for i in CID :
        analysisFolder (i)
        cenMass (i)
        cenMassSym (i)
        KMEANScluster (i)
        confWritePocket (i)
        pocketDocking (i)
        resultOutput (i)