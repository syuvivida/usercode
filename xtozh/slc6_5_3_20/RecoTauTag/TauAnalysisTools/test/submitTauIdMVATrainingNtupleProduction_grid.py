#!/usr/bin/env python

import TauAnalysis.Configuration.tools.eos as eos

import os
import shlex
import string
import subprocess
import time

samples = {
    'ZplusJets_madgraph' : {
        'datasetpath'                        : '/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
        'dbs_url'                            : 'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet',
        'events_per_job'                     : 20000,
        'total_number_of_events'             : 5000000,
        'type'                               : 'SignalMC'
    },
    'PPmuXptGt20Mu15' : {
        'datasetpath'                        : '/QCD_Pt_20_MuEnrichedPt_15_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
        'dbs_url'                            : 'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet',
        'events_per_job'                     : 20000,
        'total_number_of_events'             : 5000000,
        'type'                               : 'BackgroundMC'
    },
    'QCDmuEnrichedPt50to80' : {
        'datasetpath'                        : '/QCD_Pt-50to80_MuEnrichedPt5_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
        'dbs_url'                            : 'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet',
        'events_per_job'                     : 20000,
        'total_number_of_events'             : 2000000,
        'type'                               : 'BackgroundMC'
    },
    'QCDmuEnrichedPt80to120' : {
        'datasetpath'                        : '/QCD_Pt-80to120_MuEnrichedPt5_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
        'dbs_url'                            : 'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet',
        'events_per_job'                     : 20000,
        'total_number_of_events'             : 2000000,
        'type'                               : 'BackgroundMC'
    },
    'QCDmuEnrichedPt120to170' : {
        'datasetpath'                        : '/QCD_Pt-120to170_MuEnrichedPt5_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
        'dbs_url'                            : 'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet',
        'events_per_job'                     : 20000,
        'total_number_of_events'             : 2000000,
        'type'                               : 'BackgroundMC'
    },    
    'QCDmuEnrichedPt170to300' : {
        'datasetpath'                        : '/QCD_Pt-170to300_MuEnrichedPt5_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
        'dbs_url'                            : 'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet',
        'events_per_job'                     : 20000,
        'total_number_of_events'             : 2000000,
        'type'                               : 'BackgroundMC'
    },
    'QCDmuEnrichedPt300to470' : {
        'datasetpath'                        : '/QCD_Pt-300to470_MuEnrichedPt5_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
        'dbs_url'                            : 'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet',
        'events_per_job'                     : 20000,
        'total_number_of_events'             : 2000000,
        'type'                               : 'BackgroundMC'
    },
    'QCDmuEnrichedPt470to600' : {
        'datasetpath'                        : '/QCD_Pt-470to600_MuEnrichedPt5_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
        'dbs_url'                            : 'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet',
        'events_per_job'                     : 20000,
        'total_number_of_events'             : 3000000,
        'type'                               : 'BackgroundMC'
    },
    'QCDmuEnrichedPt600to800' : {
        'datasetpath'                        : '/QCD_Pt-600to800_MuEnrichedPt5_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
        'dbs_url'                            : 'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet',
        'events_per_job'                     : 20000,
        'total_number_of_events'             : 4000000,
        'type'                               : 'BackgroundMC'
    },
    'QCDmuEnrichedPt800to1000' : {
        'datasetpath'                        : '/QCD_Pt-800to1000_MuEnrichedPt5_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
        'dbs_url'                            : 'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet',
        'events_per_job'                     : 20000,
        'total_number_of_events'             : 3000000,
        'type'                               : 'BackgroundMC'
    },
    'QCDmuEnrichedPtGt1000' : {
        'datasetpath'                        : '/QCD_Pt-1000_MuEnrichedPt5_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
        'dbs_url'                            : 'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet',
        'events_per_job'                     : 20000,
        'total_number_of_events'             : 2000000,
        'type'                               : 'BackgroundMC'
    },
    'WplusJets_madgraph' : {
        'datasetpath'                        : '/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v2/AODSIM',
        'dbs_url'                            : 'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet',
        'events_per_job'                     : 20000,
        'total_number_of_events'             : 5000000,
        'type'                               : 'BackgroundMC'
    },
    'QCDjetsFlatPt15to3000' : {
        'datasetpath'                        : '/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
        'dbs_url'                            : 'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet',
        'events_per_job'                     : 20000,
        'total_number_of_events'             : 9991674,
        'type'                               : 'BackgroundMC'
    },
    'QCDjetsPt50to80' : {
        'datasetpath'                        : '/QCD_Pt-50to80_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v2/AODSIM',
        'dbs_url'                            : 'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet',
        'events_per_job'                     : 20000,
        'total_number_of_events'             : 2000000,
        'type'                               : 'BackgroundMC'
    },
    'QCDjetsPt80to120' : {
        'datasetpath'                        : '/QCD_Pt-80to120_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v3/AODSIM',
        'dbs_url'                            : 'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet',
        'events_per_job'                     : 20000,
        'total_number_of_events'             : 2000000,
        'type'                               : 'BackgroundMC'
    },
    'QCDjetsPt120to170' : {
        'datasetpath'                        : '/QCD_Pt-120to170_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v3/AODSIM',
        'dbs_url'                            : 'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet',
        'events_per_job'                     : 20000,
        'total_number_of_events'             : 2000000,
        'type'                               : 'BackgroundMC'
    },        
    'QCDjetsPt170to300' : {
        'datasetpath'                        : '/QCD_Pt-170to300_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v2/AODSIM',
        'dbs_url'                            : 'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet',
        'events_per_job'                     : 20000,
        'total_number_of_events'             : 2000000,
        'type'                               : 'BackgroundMC'
    },
    'QCDjetsPt300to470' : {
        'datasetpath'                        : '/QCD_Pt-300to470_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v2/AODSIM',
        'dbs_url'                            : 'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet',
        'events_per_job'                     : 20000,
        'total_number_of_events'             : 2000000,
        'type'                               : 'BackgroundMC'
    },
    'QCDjetsPt470to600' : {
        'datasetpath'                        : '/QCD_Pt-470to600_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v2/AODSIM',
        'dbs_url'                            : 'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet',
        'events_per_job'                     : 20000,
        'total_number_of_events'             : 3000000,
        'type'                               : 'BackgroundMC'
    },
    'QCDjetsPt600to800' : {
        'datasetpath'                        : '/QCD_Pt-600to800_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v2/AODSIM',
        'dbs_url'                            : 'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet',
        'events_per_job'                     : 20000,
        'total_number_of_events'             : 4000000,
        'type'                               : 'BackgroundMC'
    },
    'QCDjetsPt800to1000' : {
        'datasetpath'                        : '/QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v2/AODSIM',
        'dbs_url'                            : 'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet',
        'events_per_job'                     : 20000,
        'total_number_of_events'             : 3000000,
        'type'                               : 'BackgroundMC'
    },
    'QCDjetsPt1000to1400' : {
        'datasetpath'                        : '/QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
        'dbs_url'                            : 'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet',
        'events_per_job'                     : 20000,
        'total_number_of_events'             : 2000000,
        'type'                               : 'BackgroundMC'
    },
    'QCDjetsPt1400to1800' : {
        'datasetpath'                        : '/QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
        'dbs_url'                            : 'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet',
        'events_per_job'                     : 20000,
        'total_number_of_events'             : 2000000,
        'type'                               : 'BackgroundMC'
    },
    'QCDmuEnrichedPtGt1800' : {
        'datasetpath'                        : '/QCD_Pt-1800_TuneZ2star_8TeV_pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM',
        'dbs_url'                            : 'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet',
        'events_per_job'                     : 20000,  
        'total_number_of_events'             : 2000000,
        'type'                               : 'BackgroundMC'
    }
}
smHiggsMassPoints = [ 80, 90, 100, 110, 120, 130, 140 ]
for massPoint in smHiggsMassPoints:
    ggSampleName = "ggHiggs%1.0ftoTauTau" % massPoint
    samples[ggSampleName] = {
        'datasetpath'                        : '/GluGluToHToTauTau_M-%1.0f_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM' % massPoint,
        'dbs_url'                            : 'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet',
        'events_per_job'                     : 20000,
        'total_number_of_events'             : 500000,
        'type'                               : 'SignalMC'
    }
    vbfSampleName = "vbfHiggs%1.0ftoTauTau" % massPoint
    samples[vbfSampleName] = {
        'datasetpath'                        : '/VBF_HToTauTau_M-%1.0f_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM' % massPoint,
        'dbs_url'                            : 'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet',
        'events_per_job'                     : 20000,
        'total_number_of_events'             : 500000,
        'type'                               : 'SignalMC'
    }
mssmHiggsMassPoints = [ 160, 180, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000 ]
for massPoint in mssmHiggsMassPoints:
    ggSampleName = "ggA%1.0ftoTauTau" % massPoint
    samples[ggSampleName] = {
        'datasetpath'                        : '/SUSYGluGluToHToTauTau_M-%1.0f_8TeV-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM' % massPoint,
        'dbs_url'                            : 'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet',
        'events_per_job'                     : 20000,
        'total_number_of_events'             : 200000,
        'type'                               : 'SignalMC'
    }
    bbSampleName = "bbA%1.0ftoTauTau" % massPoint
    samples[bbSampleName] = {
        'datasetpath'                        : '/SUSYBBHToTauTau_M-%1.0f_8TeV-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM' % massPoint,
        'dbs_url'                            : 'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet',
        'events_per_job'                     : 20000,
        'total_number_of_events'             : 200000,
        'type'                               : 'SignalMC'
    }
ZprimeMassPoints = [ 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500 ]
for massPoint in ZprimeMassPoints:
    sampleName = "Zprime%1.0ftoTauTau" % massPoint
    samples[sampleName] = {
        'datasetpath'                        : '/ZprimeSSMToTauTau_M-%1.0f_TuneZ2star_8TeV-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM' % massPoint,
        'dbs_url'                            : 'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet',
        'events_per_job'                     : 20000,
        'total_number_of_events'             : -1,
        'type'                               : 'SignalMC'
    }
WprimeMassPoints = [ 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000, 3200, 3500, 4000 ]
for massPoint in WprimeMassPoints:
    sampleName = "Wprime%1.0ftoTauNu" % massPoint
    samples[sampleName] = {
        'datasetpath'                        : '/WprimeToTauNu_M-%1.0f_TuneZ2star_8TeV-pythia6-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM' % massPoint,
        'dbs_url'                            : 'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet',
        'events_per_job'                     : 20000,
        'total_number_of_events'             : -1,
        'type'                               : 'SignalMC'
    }

version = "tauId_v1_14"

submitJobFraction = 1.00

crab_template_mc = string.Template('''
[CRAB]
jobtype = cmssw
scheduler = remoteGlidein
use_server = 0

[CMSSW]
datasetpath = $datasetpath
dbs_url = $dbs_url
pset = $pset
output_file = tauIdMVATrainingNtuple.root
total_number_of_events = $total_number_of_events
events_per_job = $events_per_job

[USER]
ui_working_dir = $ui_working_dir
return_data = 0
copy_data = 1
publish_data = 0
storage_element = T2_CH_CERN
user_remote_dir = $user_remote_dir
debug_wrapper = 1

[GRID]
##SE_white_list = T2_DE_DESY
SE_black_list = T2_US_Nebraska,T2_KR_KNU,T2_IT_Legnaro
''')

crab_template_data = string.Template('''
[CRAB]
jobtype = cmssw
scheduler = remoteGlidein
use_server = 0
 
[CMSSW]
datasetpath = $datasetpath
dbs_url = $dbs_url
pset = $pset
output_file = $output_file
lumi_mask = /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt
total_number_of_lumis = -1
lumis_per_job = $lumis_per_job
#runselection = 190450-190790

[USER]
ui_working_dir = $ui_working_dir
return_data = 0
copy_data = 1
publish_data = 0
storage_element = T2_CH_CERN
user_remote_dir = $user_remote_dir
debug_wrapper = 1

[GRID]
##SE_white_list = T2_DE_DESY
SE_black_list = T2_US_Nebraska,T2_KR_KNU,T2_IT_Legnaro
''')

configFile = "produceTauIdMVATrainingNtuple_cfg.py"

currentDirectory    = os.getcwd()
submissionDirectory = os.path.join(currentDirectory, "crab")

executable_crab = 'crab'
#executable_crab = 'crab -GRID.dont_check_proxy 1' # NOTE: requires to execute 'voms-proxy-init -voms cms -valid 72:0' prior to running submitAntiMuonDiscrMVATrainingNtupleProduction_grid.py

def getStringRep_bool(flag):
    retVal = None
    if flag:
        retVal = "True"
    else:
        retVal = "False"
    return retVal

def runCommand(commandLine):
    print(commandLine)
    subprocess.call(commandLine, shell = True)

def createFilePath(filePath):
    try:
        eos.lsl(filePath)
    except IOError:
        print "filePath = %s does not yet exist, creating it." % filePath
        eos.mkdir(filePath)
        time.sleep(3)
    eos.chmod(filePath, 777)

crabCommands_create_and_submit = []
crabCommands_publish           = []

for sampleName, sampleOption in samples.items():
    
    # create config file for cmsRun
    cfgFileName_original = configFile
    cfgFile_original = open(cfgFileName_original, "r")
    cfg_original = cfgFile_original.read()
    cfgFile_original.close()

    cfg_modified = cfg_original.replace("#__", "")
    cfg_modified = cfg_modified.replace("#type#", "'%s'" % sampleOption['type'])

    cfgFileName_modified = os.path.join(submissionDirectory, cfgFileName_original.replace("_cfg.py", "_%s_%s_cfg.py" % (sampleName, version)))
    cfgFile_modified = open(cfgFileName_modified, "w")
    cfgFile_modified.write(cfg_modified)
    cfgFile_modified.close()

    output_files = [ "tauIdMVATrainingNtuple.root" ]
        
    # create crab config file
    crabOptions = None
    crab_template = None    
    if sampleOption['type'] == "SignalMC" or sampleOption['type'] == "BackgroundMC":
        total_number_of_events = None
        if submitJobFraction < 1.0 and sampleOption['total_number_of_events'] != -1:
            print "submitting fraction = %1.2f of event statistics for sample = %s" % (submitJobFraction, sampleName)
            total_number_of_events = int(submitJobFraction*sampleOption['total_number_of_events'])
        else:
            total_number_of_events = sampleOption['total_number_of_events']
        crabOptions = {
            'datasetpath'            : sampleOption['datasetpath'],
            'dbs_url'                : sampleOption['dbs_url'],
            'total_number_of_events' : total_number_of_events,
            'events_per_job'         : sampleOption['events_per_job'],
            'pset'                   : cfgFileName_modified,
            'output_file'            : ",".join(output_files),
            'ui_working_dir'         : os.path.join(submissionDirectory, "crabdir_%s_%s" % (sampleName, version)),
            'user_remote_dir'        : "CMSSW_5_3_x/Ntuples/tauIdMVATraining/%s/%s" % (version, sampleName)
        }
        crab_template = crab_template_mc
    elif sampleOption['type'] == "Data":
        crabOptions = {
            'datasetpath'            : sampleOption['datasetpath'],
            'dbs_url'                : sampleOption['dbs_url'],
            'lumis_per_job'          : sampleOption['lumis_per_job'],
            'pset'                   : cfgFileName_modified,
            'output_file'            : ",".join(output_files),
            'ui_working_dir'         : os.path.join(submissionDirectory, "crabdir_%s_%s" % (sampleName, version)),
            'user_remote_dir'        : "CMSSW_5_3_x/Ntuples/tauIdMVATraining/%s/%s" % (version, sampleName)
        }
        crab_template = crab_template_data
    else:
        raise ValueError("Invalid sample type = %s !!" % sampleOption['type'])
    crabFileName = "crab_tauIdMVATrainingNtupleProduction_%s_%s.cfg" % (sampleName, version)
    crabFileName_full = os.path.join(submissionDirectory, crabFileName)
    crabFile = open(crabFileName_full, 'w')
    crabConfig = crab_template.substitute(crabOptions)
    crabFile.write(crabConfig)
    crabFile.close()

    # create output directory
    # (in principle crab will do this, but sometimes fails with 'Permission denied' error, causing all jobs to fail with error code 60307)
    createFilePath("/store/user/veelken/CMSSW_5_3_x/Ntuples/tauIdMVATraining/%s" % version)
    createFilePath("/store/user/veelken/CMSSW_5_3_x/Ntuples/tauIdMVATraining/%s/%s" % (version, sampleName))

    # keep track of commands necessary to create, submit and publish crab jobs
    crabCommands_create_and_submit.append('%s -create -cfg %s' % (executable_crab, crabFileName_full))
    if 'events_per_job' in sampleOption.keys(): # MC
        events_total = None
        if sampleOption['total_number_of_events'] == -1:
            events_total = 10000000
        else:
            events_total = sampleOption['total_number_of_events']
        if (events_total / sampleOption['events_per_job']) < 450: # CV: add 10% safety margin to avoid jobs not getting submitted at all in case crab decides to create more than 500 jobs
            crabCommands_create_and_submit.append('%s -submit -c %s' % (executable_crab, crabOptions['ui_working_dir']))
        else:
            numJobs = (events_total / sampleOption['events_per_job'])
            if (events_total % sampleOption['events_per_job']) != 0:
                numJobs = numJobs + 1
            numJobs_per_submitCall = 500
            numSubmitCalls = (numJobs / numJobs_per_submitCall)
            if (numJobs % numJobs_per_submitCall) != 0:
                numSubmitCalls = numSubmitCalls + 1
            for submitIdx in range(numSubmitCalls):
                jobId_first = submitIdx*500 + 1
                jobId_last  = (submitIdx + 1)*500
                if jobId_last > numJobs:
                    jobId_last = numJobs
                crabCommands_create_and_submit.append('echo "pausing for 10 seconds before submitting next batch of jobs..."')
                crabCommands_create_and_submit.append('sleep 10')
                crabCommands_create_and_submit.append('%s -submit %i-%i -c %s' % (executable_crab, jobId_first, jobId_last, crabOptions['ui_working_dir']))
    else: # Data
        crabCommands_create_and_submit.append('%s -submit -c %s' % (executable_crab, crabOptions['ui_working_dir']))
    
shellFileName_create_and_submit = "tauIdMVATrainingNtupleProduction_crab_create_and_submit.sh"
shellFile_create_and_submit = open(shellFileName_create_and_submit, "w")
for crabCommand in crabCommands_create_and_submit:
    shellFile_create_and_submit.write("%s\n" % crabCommand)
shellFile_create_and_submit.close()

print("Finished building config files. Now execute 'source %s' to create & submit crab jobs." % shellFileName_create_and_submit)

