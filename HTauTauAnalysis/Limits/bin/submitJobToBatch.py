#!/usr/bin/env python

import commands
import re
import os

import sys
sys.path.append('./')

###########################################
###########################################

def treeSkim( sample, xsection, runInSeries=True):

    os.system( 'mkdir batch/' )

    stream = "MuTau"
    #if(re.search("MuTau",sample)!=None):
    #    stream = "MuTau"
    #else:
    #    stream = "ElecTau"
    print "Stream ", stream
    
    if runInSeries:
         print "Running in series via the command ..."
         runInSeries   = 'treeSkimmer'+stream+'_MVA '+sample+' '+str(xsection)
         print runInSeries
         os.system(runInSeries)
         mv     = 'mv nTuple'+sample+'_Open_'+stream+'Stream_*.root batch/' 
         hadd   = 'hadd -f batch/nTuple'+sample+'_Open_'+stream+'Stream.root batch/nTuple'+sample+'_Open_'+stream+'Stream_*.root'
         remove = 'rm batch/nTuple'+sample+'_Open_'+stream+'Stream_*.root'
         print 'Now doing hadd...'
         os.system(mv)
         os.system(hadd)
         os.system(remove)

    else:
        print "Creating the shell file for the batch..."
        f = open('batch/job'+'_'+sample+'.sh','w')    
        f.write('#!/bin/sh\n\n')
        f.write('export WORKINGDIR="/home/llr/cms/lbianchini/CMSSW_4_2_8_patch7_reload/src/Bianchi/Limits/bin/"\n')
        f.write('source /opt/exp_soft/cms/cmsset_default.sh\n')
        f.write('cd $WORKINGDIR\n')
        f.write('eval `scram runtime -sh`\n')
        f.write('source /opt/exp_soft/cms/cms_ui_env_crab.sh\n')
        f.write('treeSkimmer'+stream+'_MVA '+sample+' '+str(xsection)+' \n')
        f.write('mv nTuple'+sample+'_Open_'+stream+'Stream_*.root batch/\n')
        f.write('hadd -f batch/nTuple'+sample+'_Open_'+stream+'Stream.root batch/nTuple'+sample+'_Open_'+stream+'Stream_*.root\n')
        f.write('rm batch/nTuple'+sample+'_Open_'+stream+'Stream_*.root\n')
        f.close()
        #os.system('source /opt/exp_soft/cms/t3/t3setup')
        print "Submitting the job to the cms queue via the command ..."
        submitToQueue = '/opt/exp_soft/cms/t3/t3submit -q llr -k eo -N '+sample+' /home/llr/cms/lbianchini/CMSSW_4_2_8_patch7_reload/src/Bianchi/Limits/bin/batch/job'+'_'+sample+'.sh > batch/'+sample+'.txt'
        print submitToQueue
        os.system(submitToQueue)
  
###########################################
###########################################


#treeSkim("GGFH125",      6.37e-02*15.31 * 1.0 * 0.0589651)
treeSkim("VBFH125",      6.37e-02*1.211 * 1.0 * 0.0890235)

#treeSkim("GGFH125",      6.37e-02*15.31 * 1.0 * 0.037848)
#treeSkim("VBFH125",      6.37e-02*1.211 * 1.0 * 0.061610) 


#treeSkim("Run2011A-MuTau-05AugReReco-v2-Embedded-t3_run",     0)   

#treeSkim("Run2011-MuTau-All_run",              0)                            
#treeSkim("Run2011-MuTau-Embedded-All_run",     0)                           

#treeSkim("DYJets-MuTau-50-madgraph-PUS6_run",  3048           * 0.009631    * 0.704144) 
#treeSkim("TTJets-MuTau-madgraph-PUS6_run",     157.5          * 0.020998    * 0.875478)  
#treeSkim("WJets-MuTau-madgraph-PUS6_run",      31314.0        * 0.001261    * 0.645690)  
#treeSkim("W3Jets-MuTau-madgraph-PUS6_run",     304.0          * 1.0         * 0.13987)
#treeSkim("WZ-MuTau-pythia-PUS6_run",           18.2           * 0.0068962   * 0.8181)        
#treeSkim("ZZ-MuTau-pythia-PUS6_run",            5.9           * 0.0060357   * 0.814714)    
#treeSkim("WW-MuTau-pythia-PUS6_run",           43.0           * 1.0         * 0.075776)       

#treeSkim("GGFH110-MuTau-powheg-PUS6_run",      8.02e-02*19.84 * 1.0         * 0.074268)    
#treeSkim("GGFH115-MuTau-powheg-PUS6_run",      7.65e-02*18.13 * 1.0         * 0.080647)  
#treeSkim("GGFH120-MuTau-powheg-PUS6_run",      7.10e-02*16.63 * 1.0         * 0.084353)  
#treeSkim("GGFH125-MuTau-powheg-PUS6_run",      6.37e-02*15.31 * 1.0         * 0.088036) 
#treeSkim("GGFH130-MuTau-powheg-PUS6_run",      5.48e-02*14.12 * 0.05479     * 0.892072)
#treeSkim("GGFH135-MuTau-powheg-PUS6_run",      4.52e-02*13.08 * 1.0         * 0.097084)    
#treeSkim("GGFH140-MuTau-powheg-PUS6_run",      3.54e-02*12.13 * 0.05808     * 0.889078)  
#treeSkim("GGFH145-MuTau-powheg-PUS6_run",      2.61e-02*11.27 * 1.0         * 0.104852) 

#treeSkim("VBFH110-MuTau-powheg-PUS6_run",      8.02e-02*1.398 * 1.0         * 0.120909)
#treeSkim("VBFH115-MuTau-powheg-PUS6_run",      7.65e-02*1.332 * 0.05825     * 0.872899)
#treeSkim("VBFH120-MuTau-powheg-PUS6_run",      7.10e-02*1.269 * 1.0         * 0.128796)   
#treeSkim("VBFH125-MuTau-powheg-PUS6_run",      6.37e-02*1.211 * 0.06337     * 0.876146)   
#treeSkim("VBFH130-MuTau-powheg-PUS6_run",      5.48e-02*1.154 * 0.06515     * 0.873742)
#treeSkim("VBFH135-MuTau-powheg-PUS6_run",      4.52e-02*1.100 * 0.06740     * 0.879899)
#treeSkim("VBFH140-MuTau-powheg-PUS6_run",      3.54e-02*1.052 * 1.0         * 0.140405)
#treeSkim("VBFH145-MuTau-powheg-PUS6_run",      2.61e-02*1.004 * 0.07349     * 0.876581)

#treeSkim("VH110-MuTau-pythia-PUS6_run",        8.02e-02*(0.8754+0.4721+0.1257  ) * 0.088174 * 0.905001) 
#treeSkim("VH115-MuTau-pythia-PUS6_run",        7.65e-02*(0.7546+0.4107+0.1106  ) * 0.092032 * 0.904303)
#treeSkim("VH120-MuTau-pythia-PUS6_run",        7.10e-02*(0.6561+0.3598+0.09756 ) * 1.0      * 0.185308)   
#treeSkim("VH125-MuTau-pythia-PUS6_run",        6.37e-02*(0.5729+0.3158+0.08634 ) * 0.09885  * 0.905322)
#treeSkim("VH130-MuTau-pythia-PUS6_run",        5.48e-02*(0.4942+0.2778+0.07658 ) * 1.0      * 0.19413)
#treeSkim("VH135-MuTau-pythia-PUS6_run",        4.52e-02*(0.4390+0.2453+0.06810 ) * 0.104662 * 0.906734)
#treeSkim("VH140-MuTau-pythia-PUS6_run",        3.54e-02*(0.3857+0.2172+0.06072 ) * 0.108278 * 0.905325)  
#treeSkim("VH145-MuTau-pythia-PUS6_run",        2.61e-02*(0.3406+0.1930+0.05435 ) * 1.0      * 0.209387)  

#treeSkim("VBFH100-MuTau-powheg-PUS6_run",      8.36e-02*1.546 * 0.04953     * 0.809326) 
#treeSkim("GGFH100-MuTau-powheg-PUS6_run",      8.36e-02*24.02 * 0.03816     * 0.818372) 
#treeSkim("VBFH105-MuTau-powheg-PUS6_run",      8.25e-02*1.472 * 0.05278     * 0.804249) 
#treeSkim("GGFH105-MuTau-powheg-PUS6_run",      8.25e-02*21.78 * 0.04137     * 0.814457) 
#treeSkim("VBFH160-MuTau-powheg-PUS6_run",      5.32e-04*0.8787* 1.0         * 0.137951)  
#treeSkim("GGFH160-MuTau-powheg-PUS6_run",      5.32e-04*9.080 * 1.0         * 0.104164)  
#treeSkim("VH100-MuTau-pythia-PUS6_run",        2.61e-02*(1.186+ 0.6313+0.1638  ) * 0.082128 * 0.834434)  
#treeSkim("VH160-MuTau-pythia-PUS6_run",        5.32e-04*(0.2291+0.1334+0.03942 ) * 1.0      * 0.203857)  



###########################################
###########################################


#treeSkim("Run2011-ElecTau-All_run",              0)                            
#treeSkim("Run2011-ElecTau-Embedded-All_run",     0)                           

#treeSkim("DYJets-ElecTau-50-madgraph-PUS6_run",  3048           * 0.0537207    * 0.864464) 
#treeSkim("TTJets-ElecTau-madgraph-PUS6_run",     157.5          * 0.0149329    * 0.796639)  
#treeSkim("WJets-ElecTau-madgraph-PUS6_run",      31314.0        * 0.0011910    * 0.624904)  
##treeSkim("W3Jets-ElecTau-madgraph-PUS6_run",     304.0          * 1.0         * 0.1257)
#treeSkim("WZ-ElecTau-pythia-PUS6_run",           18.2           * 0.0111572   * 0.85868)        
#treeSkim("ZZ-ElecTau-pythia-PUS6_run",            5.9           * 0.0160517   * 0.88701)    
#treeSkim("WW-ElecTau-pythia-PUS6_run",           43.0           * 0.0092734   * 0.81103)       

#treeSkim("GGFH110-ElecTau-powheg-PUS6_run",      8.02e-02*19.84 * 0.044102646 * 0.839456)    
#treeSkim("GGFH115-ElecTau-powheg-PUS6_run",      7.65e-02*18.13 * 0.04634647  * 0.846653)  
#treeSkim("GGFH120-ElecTau-powheg-PUS6_run",      7.10e-02*16.63 * 1.0         * 0.083277)  
#treeSkim("GGFH125-ElecTau-powheg-PUS6_run",      6.37e-02*15.31 * 1.0         * 0.087446) 
#treeSkim("GGFH130-ElecTau-powheg-PUS6_run",      5.48e-02*14.12 * 1.0         * 0.091171)
#treeSkim("GGFH135-ElecTau-powheg-PUS6_run",      4.52e-02*13.08 * 0.0573246   * 0.849161)    
#treeSkim("GGFH140-ElecTau-powheg-PUS6_run",      3.54e-02*12.13 * 1.0         * 0.101619)  
#treeSkim("GGFH145-ElecTau-powheg-PUS6_run",      2.61e-02*11.27 * 1.0         * 0.104718) 

#treeSkim("VBFH110-ElecTau-powheg-PUS6_run",      8.02e-02*1.398 * 0.05580703  * 0.832940)
#treeSkim("VBFH115-ElecTau-powheg-PUS6_run",      7.65e-02*1.332 * 0.05869278  * 0.839389)
#treeSkim("VBFH120-ElecTau-powheg-PUS6_run",      7.10e-02*1.269 * 1.0         * 0.130759)   
#treeSkim("VBFH125-ElecTau-powheg-PUS6_run",      6.37e-02*1.211 * 1.0         * 0.133172)   
#treeSkim("VBFH130-ElecTau-powheg-PUS6_run",      5.48e-02*1.154 * 1.0         * 0.137893)
#treeSkim("VBFH135-ElecTau-powheg-PUS6_run",      4.52e-02*1.100 * 1.0         * 0.142746)
#treeSkim("VBFH140-ElecTau-powheg-PUS6_run",      3.54e-02*1.052 * 1.0         * 0.144337)
#treeSkim("VBFH145-ElecTau-powheg-PUS6_run",      2.61e-02*1.004 * 1.0         * 0.148619)

#treeSkim("VH110-ElecTau-pythia-PUS6_run",        8.02e-02*(0.8754+0.4721+0.1257  ) * 0.08438633     * 0.857899) 
#treeSkim("VH115-ElecTau-pythia-PUS6_run",        7.65e-02*(0.7546+0.4107+0.1106  ) * 1.0            * 0.168881)
#treeSkim("VH120-ElecTau-pythia-PUS6_run",        7.10e-02*(0.6561+0.3598+0.09756 ) * 1.0            * 0.1751)   
#treeSkim("VH125-ElecTau-pythia-PUS6_run",        6.37e-02*(0.5729+0.3158+0.08634 ) * 1.0            * 0.1793736)
#treeSkim("VH130-ElecTau-pythia-PUS6_run",        5.48e-02*(0.4942+0.2778+0.07658 ) * 0.10051        * 0.86869811)
#treeSkim("VH135-ElecTau-pythia-PUS6_run",        4.52e-02*(0.4390+0.2453+0.06810 ) * 0.101333       * 0.87368421)
#treeSkim("VH140-ElecTau-pythia-PUS6_run",        3.54e-02*(0.3857+0.2172+0.06072 ) * 1.0            * 0.19516810)  
#treeSkim("VH145-ElecTau-pythia-PUS6_run",        2.61e-02*(0.3406+0.1930+0.05435 ) * 0.10708        * 0.8707802)  
