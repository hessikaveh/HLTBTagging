
#Define the datasets the following:
#list with
#     0th element: name
#     1st element: tuple containing primary and secondary DAS dataset name
#     2nd element: 0 if Data, 1 if MC
Data = [
#    ["HLT_Ntuple_BTagging_DiLepton_v10",
#     ("/MuonEG/Run2017B-17Nov2017-v1/MINIAOD","/MuonEG/Run2017B-v1/RAW"),
#     "_RunB",
#     True],
    ["HLT_Ntuple_BTagging_DiLepton_v10",
     ("/MuonEG/Run2017C-17Nov2017-v1/MINIAOD","/MuonEG/Run2017C-v1/RAW"),
     "_RunC",
     True],
#    ["HLT_Ntuple_BTagging_DiLepton_v10",
#     ("/MuonEG/Run2017D-17Nov2017-v1/MINIAOD","/MuonEG/Run2017D-v1/RAW"),
#     "_RunD",
#     True],
#    ["HLT_Ntuple_BTagging_DiLepton_v10",
#     ("/MuonEG/Run2017E-17Nov2017-v1/MINIAOD","/MuonEG/Run2017E-v1/RAW"),
#     "_RunE",
#     True],
#    ["HLT_Ntuple_BTagging_DiLepton_v10",
#     ("/MuonEG/Run2017F-17Nov2017-v1/MINIAOD","/MuonEG/Run2017F-v1/RAW"),
#     "_RunF",
#     True]
]
MC = [
    ["HLT_Ntuple_BTagging_DiLepton_v10",
     #("/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer17DRStdmix-NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v2/AODSIM",
     ("/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAOD-TSG_94X_mc2017_realistic_v11-v1/MINIAODSIM",
      "/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIFall17DRPremix-TSG_94X_mc2017_realistic_v11-v1/GEN-SIM-RAW"),
     "",
     False],
#    ["HLT_Ntuple_BTagging_DiLepton_v10",
#     #("/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer17DRStdmix-NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v2/AODSIM",
#     ("/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAOD-TSG_94X_mc2017_realistic_v11-v1/MINIAODSIM",
#      "/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17DRPremix-TSG_94X_mc2017_realistic_v11-v1/GEN-SIM-RAW"),
#     "",
#     False],
#    ["HLT_Ntuple_BTagging_DiLepton_v10",
#     ("/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
#      "/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v2/GEN-SIM-RAW"),
#     "",
#     False]
]


#datasets = MC
datasets = Data
#datasets = [MC[1],MC[2]]
#datasets = [Data[0]]
print datasets
raw_input("press ret to continue")
prefix = "AddTrks2p2"


if __name__ == '__main__':
    from CRABAPI.RawCommand import crabCommand
    from CRABClient.UserUtilities import config
    config = config()
    
    for dataset in datasets:
        name = dataset[0]
        config.section_("General")
        config.General.workArea = 'crab_' + name + prefix 
        config.General.transferLogs=True
#       config.General.requestName = name+"_"+dataset.replace('/',"_")
        config.General.requestName = name + "_V3_"  + prefix + "_" + dataset[1][0].split('/')[1].split("-")[0] + dataset[2]

        print "Requestname: ", name + "_V3_" + prefix + "_" + dataset[1][0].split('/')[1].split("-")[0] + dataset[2]
        #raw_input("press ret")
        
        config.section_("JobType")
#        config.JobType.numCores = 4
        config.JobType.numCores = 4
        config.JobType.maxMemoryMB = 5000
        config.JobType.maxJobRuntimeMin = 2000
        config.JobType.sendPythonFolder = True
        config.JobType.pluginName = 'Analysis'
        config.JobType.psetName = 'crab_fake_pset.py'
        config.JobType.scriptExe = 'crab_script_phaseI.sh'
        import os
        os.system("tar czf python.tar.gz --dereference --directory $CMSSW_BASE python")
        os.system("voms-proxy-info -path | xargs -i  cp {}  .")
        config.JobType.inputFiles = [
            'hlt_dump_phase1.py',
            'hlt_dump_mc_phase1.py',
            'fwlite_config_phaseI.py',
            'script_phaseI.py',
            'utils.py',
            'python.tar.gz',
        ]
        
        config.section_("Data")
        config.Data.inputDBS = 'global'
        #Data
        if dataset[3]:
            config.Data.splitting = 'LumiBased'
            config.Data.unitsPerJob = 70##FIXME: use 20
        #MC
        if not dataset[3]:
            config.Data.splitting = 'LumiBased'
            config.Data.unitsPerJob = 35
            ##FIXME: use 20
        

        config.Data.totalUnits = -1 #10*config.Data.unitsPerJob #FIXME: use -1
        #config.Data.outLFNDirBase = '/store/user/koschwei/onlineBTV/' + name + prefix
        config.Data.outLFNDirBase = '/store/user/hkaveh/' + name + prefix
        config.Data.publication = True
        if dataset[3]:
            print "Using JSON"
            config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'
            #config.Data.lumiMask = '/afs/cern.ch/work/k/koschwei/public/test/CMSSW_9_2_12_patch1/src/HLTBTagging/nTuples/PU28to63_Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
        config.Data.inputDataset = dataset[1][0]
        config.Data.secondaryInputDataset = dataset[1][1]
#       config.Data.publishDataName = config.General.requestName
        config.Data.outputDatasetTag = name
        config.Data.allowNonValidInputDataset = True
        #config.Site.blacklist = ['T0_*','T3_*','T2_AT_*','T2_CH_*','T2_BE_*','T2_BR_*','T2_DE_RWTH','T2_ES_*','T2_EE_*','T2_FR_*','T2_FI_*','T2_GR_*','T2_HU_*','T2_IN_*','T2_IT_*','T2_KR_*','T2_MY_*','T2_PK_*','T2_PL_*','T2_PT_*','T2_RU_*','T2_TH_*','T2_TR_*','T2_TW_*','T2_UA_*','T2_UK_*','T2_US_*'] # == config.Site.whitelist = ["T2_DE_DESY"]
        #config.Site.whitelist = ["T2_DE_DESY"]

        
        config.section_("Site")
        #config.Site.storageSite = "T2_CH_CSCS"
        config.Site.storageSite = 'T2_DE_DESY' 
        #config.Site.storageSite = "T3_CH_PSI"
        print "submitting ",dataset
        crabCommand('submit',config = config)
        print "DONE ",dataset
    
