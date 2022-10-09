from distutils.ccompiler import CCompiler
#from tty import CC
from chimerax.core.commands import run
import os
import numpy as np

PAATH = 'C:/Users/c/Dropbox/27. 연구 (2022)/__TMEM16A & Drug-Binding Site 연구__/pocketAnalysis'

def OOOpen () :
    global receptor, PAATH
    run (session, "close all")
    run (session, "open '%s/%s.pdb'" % (PAATH, receptor))
    run (session, "view orient; view initial")
    return

def CCColor () :
    run (session, "color steel blue")
    run (session, "color byelement")
    run (session, "turn x 90")
    run (session, "turn z -1")
    run (session, "windowsize 1200 800")
    run (session, "set bgColor white")
    run (session, "graphics silhouettes true")
    run (session, "lighting soft")
    return


def chimeraView (cid) :
    global receptor, PAATH

    #run (session, "open '%s/analysis_%s/pocketDocking__final.pdb'" % (PAATH, cid))


    # docking
    OOOpen ()
    for file in os.listdir ("%s/linux_%s/outResult" % (PAATH, cid)):
        run (session, "open '%s/linux_%s/outResult/%s'" % (PAATH, cid, file))
    CCColor ()
    run (session, "save '%s/%s_1_docking.jpeg'" % (PAATH, cid))
    run (session, "close all")

    # cenMass
    OOOpen ()
    CCColor ()
    f = open ("%s/analysis_%s/cenMass.txt" % (PAATH, cid), "r")
    for line in f :
        line = line.strip ().split ("\t")
        temp = []
        for i in line :
            temp.append (float (i))
        [X, Y, Z] = temp
        run (session, "shape sphere center %f,%f,%f radius 1 color light coral" % (X, Y, Z))
    f.close ()
    run (session, "save '%s/%s_2_cenMass.jpeg'" % (PAATH, cid))
    run (session, "close all")


    # cenMassSym
    OOOpen ()
    CCColor ()
    f = open ("%s/analysis_%s/cenMass_sym.txt" % (PAATH, cid), "r")
    for line in f :
        line = line.strip ().split ("\t")
        temp = []
        for i in line :
            temp.append (float (i))
        [X, Y, Z] = temp
        run (session, "shape sphere center %f,%f,%f radius 1 color cornflower blue" % (X, Y, Z))
    f.close ()
    run (session, "save '%s/%s_3_cenMassSym.jpeg'" % (PAATH, cid))
    run (session, "close all")


    # clustering
    OOOpen ()
    CCColor ()
    f = open ("%s/analysis_%s/cenMass_sym_cluster.txt" % (PAATH, cid), "r")
    for line in f :
        line = line.strip ().split ("\t")
        temp = []
        for i in line [2 : ] :
            temp.append (float (i))
        [X, Y, Z] = temp
        run (session, "shape sphere center %f,%f,%f radius %f color 128,255,128,0.7 mesh true" % (X, Y, Z, float (line [0]) / 30))
    f.close ()
    run (session, "save '%s/%s_4_clustering.jpeg'" % (PAATH, cid))
    run (session, "close all")


    """
    run (session, "close all")
    run (session, "open 'C:/Users/c/Dropbox/27. 연구 (2022)/__TMEM16A & Honokiol Magnolol 연구__/pocketAnalysis/%s.pdb'" % (receptor))
    run (session, "open 'C:/Users/c/Dropbox/27. 연구 (2022)/__TMEM16A & Honokiol Magnolol 연구__/pocketAnalysis/analysis_%s/pocketDocking__final.pdb'" % (cid))
    changeChains(openModels.list(modelTypes=[Molecule]), [('', 'C')])
    run (session, "combine #0,1 model #2")
    run (session, "write format pdb 2 'C:/Users/c/Dropbox/27. 연구 (2022)/__TMEM16A & Honokiol Magnolol 연구__/pocketAnalysis/analysis_%s/pocketDocking__final_combined.pdb'" % (cid))
    """

    return

receptor = "hT16A_m_dimer"
CID =   ["72300", # magnolol
        "72303", # honokiol
        "2898877", # CaCCinh-A01
        "711253" # 1PBC
        ]
for i in CID :
    chimeraView (i)
