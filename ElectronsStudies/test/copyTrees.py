#!/usr/bin/env python

import commands
import re
import os

import sys
#sys.path.append('')

###########################################
###########################################

def copyAndAdd( directory ):
    
    outSamples = re.split(r'\n',commands.getoutput("dpns-ls /dpm/in2p3.fr/home/cms/trivcat/store/user/inaranjo/"+directory))

    for sample in  outSamples:
        samplesplit =  re.split(r'/',sample)
        output = commands.getoutput("dpns-ls /dpm/in2p3.fr/home/cms/trivcat/store/user/inaranjo/"+directory+"/"+samplesplit[len(samplesplit)-1])
        
        if  re.search("DYJets-iter", samplesplit[len(samplesplit)-1])==None:
            continue


        
        outFiles = re.split(r'\n',output)
        for name in outFiles:
            splitname = re.split(r'/',name)
            print "copying file "+splitname[ len(splitname) -1 ]
            cp = "rfcp /dpm/in2p3.fr/home/cms/trivcat/store/user/inaranjo/"+directory+"/"+samplesplit[len(samplesplit)-1]+"/"+splitname[ len(splitname) -1 ]+" ."
            print cp
            os.system(cp)
        hadd = "hadd -f  MVAElecTau"+samplesplit[len(samplesplit)-1]+".root MVAElecTau_*_*_*.root"
        print hadd
        os.system(hadd)
        ###os.system("mv  treeElecTauStream_"+samplesplit[len(samplesplit)-1]+".root ../")
        remove = "rm -r MVAElecTau_*_*_*.root"
        print remove
        os.system(remove)
 
###########################################
###########################################
def removeFiles( directory ):

    outSamples = re.split(r'\n',commands.getoutput("dpns-ls /dpm/in2p3.fr/home/cms/trivcat/store/user/inaranjo/"+directory))

    for sample in  outSamples:
        samplesplit =  re.split(r'/',sample)
        output = commands.getoutput("dpns-ls /dpm/in2p3.fr/home/cms/trivcat/store/user/inaranjo/"+directory+"/"+samplesplit[len(samplesplit)-1])
        
        if  re.search(".root", samplesplit[len(samplesplit)-1])==None:
            continue

        outFiles = re.split(r'\n',output)
        for name in outFiles:
            splitname = re.split(r'/',name)
            print "removing file "+splitname[ len(splitname) -1 ]
            rm = "rfrm /dpm/in2p3.fr/home/cms/trivcat/store/user/inaranjo/"+directory+"/"+samplesplit[len(samplesplit)-1]
            print rm
            #os.system(rm)
           
        
 
###########################################
###########################################
           
## copyAndAdd( "ElecTauStreamTest" )
## removeFiles( "AntiEMVAVariables" )
copyAndAdd( "AntiEMVAVariables" )
