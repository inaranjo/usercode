#!/usr/bin/env python

import commands
import re
import os

import sys

###########################################
###########################################

def publishSkim( stream, tasks ):

    f = open(stream+'_Skim.txt', 'w')
    
    for task in tasks:
        READ = 0
        SKIM = 0
        #status = "crab -status -c "+task+"_skim/"
        status = "crab -status -c "+task
        print status
        os.system( status )
        #getoutput = "crab -getoutput -c "+task+"_skim/"
        getoutput = "crab -getoutput -c "+task
        print getoutput
        os.system( getoutput )
        #report = "crab -report -c "+task+"_skim/"
        report = "crab -report -c "+task
        print report
        output0 = commands.getoutput( report )
        reportLines = re.split(r'\n',output0)
        for line in reportLines:
            if  re.search("Total Events read:", line )!=None:
                words = re.split(r'\s', line)
                print words
                READ = float(words[3])
                
        #publish = "crab -publish -c "+task+"_skim/"
        publish = "crab -publish -c "+task
        print publish
        #os.system( publish )
        output1 = commands.getoutput( publish )
        pusblishLines = re.split(r'\n',output1)
        for line in pusblishLines:
            if  re.search("total events:", line )!=None:
                words0 = re.split(r'/',line)
                words1 = re.split(r'\s', line)
                print words1, words0
                SKIM = float(words1[3])
                datasetpath = "/"+words0[1]+"/"+words0[2]+"/"+words0[3]
                #print words
                f.write('['+task+'_run]\n')
                f.write('#skim efficiency: '+str(SKIM)+'/'+str(READ)+' = '+str(SKIM/READ)+'\n')
                f.write('CMSSW.datasetpath='+datasetpath+'\n')
                f.write('CMSSW.total_number_of_events= -1\n')
                f.write('CMSSW.events_per_job = 5000\n')
                f.write('\n')
    f.close()

###########################################
###########################################


tasksElecTauStream  = [
    #'QCD-20to30-BCtoE-ElecTau-pythia-PUS6_skim',
    #'QCD-30to80-BCtoE-ElecTau-pythia-PUS6_skim'
    'QCD-20to30-EMEnriched-ElecTau-pythia-PUS6-v4_skim'
    ]
#publishSkim( "ElecTauStream_17Feb2012-v3_EMEnriched",  tasksElecTauStream )
publishSkim( "ElecTauStream_21Feb2012_EMEnriched",  tasksElecTauStream )


## tasksElecTauStream  = [
##     'DYJets-ElecTau-50-madgraph-PUS6-v2_skim',
##     'GGFH100-ElecTau-powheg-PUS6_skim',
##     'GGFH105-ElecTau-powheg-PUS6_skim',
##     'GGFH110-ElecTau-powheg-PUS6_skim',
##     'GGFH115-ElecTau-powheg-PUS6_skim',
##     'GGFH120-ElecTau-powheg-PUS6_skim',
##     'GGFH125-ElecTau-powheg-PUS6_skim',
##     'GGFH130-ElecTau-powheg-PUS6_skim',
##     'GGFH135-ElecTau-powheg-PUS6_skim',
##     'GGFH140-ElecTau-powheg-PUS6_skim',
##     'GGFH145-ElecTau-powheg-PUS6_skim',
##     'Run2011A-ElecTau-03OctReReco_skim',
##     'Run2011A-ElecTau-05AugReReco-Embedded_skim',
##     'Run2011A-ElecTau-05AugReReco_skim',
##     'Run2011A-ElecTau-May10ReReco-Embedded_skim',
##     'Run2011A-ElecTau-May10ReReco_skim',
##     'Run2011A-ElecTau-PromptReco-v4-Embedded_skim',
##     'Run2011A-ElecTau-PromptReco-v4_skim',
##     'Run2011A-ElecTau-PromptReco-v6-Embedded_skim',
##     'Run2011A-ElecTau-PromptReco-v6-part1_skim',
##     'Run2011A-ElecTau-PromptReco-v6-part2_skim',
##     'Run2011B-ElecTau-PromptReco-v1-Embedded_skim',
##     'Run2011B-ElecTau-PromptReco-v1-part1_skim',
##     'Run2011B-ElecTau-PromptReco-v1-part2_skim',
##     'Run2011B-ElecTau-PromptReco-v1-part3_skim',
##     'Run2011B-ElecTau-PromptReco-v1-part4_skim',
##     'Run2011B-ElecTau-PromptReco-v1-part5_skim',
##     'TTJets-ElecTau-madgraph-PUS6_skim',
##     'VBFH100-ElecTau-powheg-PUS6_skim',
##     'VBFH105-ElecTau-powheg-PUS6_skim',
##     'VBFH110-ElecTau-powheg-PUS6_skim',
##     'VBFH115-ElecTau-powheg-PUS6_skim',
##     'VBFH120-ElecTau-powheg-PUS6_skim',
##     'VBFH125-ElecTau-powheg-PUS6_skim',
##     'VBFH130-ElecTau-powheg-PUS6_skim',
##     'VBFH135-ElecTau-powheg-PUS6_skim',
##     'VBFH140-ElecTau-powheg-PUS6_skim',
##     'VBFH145-ElecTau-powheg-PUS6_skim',
##     'VH100-ElecTau-pythia-PUS6_skim',
##     'VH110-ElecTau-pythia-PUS6_skim',
##     'VH115-ElecTau-pythia-PUS6_skim',
##     'VH120-ElecTau-pythia-PUS6_skim',
##     'VH125-ElecTau-pythia-PUS6_skim',
##     'VH130-ElecTau-pythia-PUS6_skim',
##     'VH135-ElecTau-pythia-PUS6_skim',
##     'VH140-ElecTau-pythia-PUS6_skim',
##     'VH145-ElecTau-pythia-PUS6_skim',
##     'W3Jets-ElecTau-madgraph-PUS6_skim',
##     'WJets-ElecTau-madgraph-PUS6_skim',
##     'WW-ElecTau-pythia-PUS6_skim',
##     'WZ-ElecTau-pythia-PUS6_skim',
##     'ZZ-ElecTau-pythia-PUS6_skim'
##     ]
## publishSkim( "ElecTauStream_08Feb2012",  tasksElecTauStream )


## tasksTest  = [
##     'VBFH115-ElecTau-powheg-PUS6_skim'
##     ]

## publishSkim( "ElecTauStream_Test",  tasksTest   )
