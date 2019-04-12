import FWCore.ParameterSet.Config as cms

process = cms.Process("FAKE")

#process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/k/koschwei/work/public/MuonEG_Run299368_PromptReco-v1_Run2017C_AOD_LS-79to90-115to129.root'),
#    secondaryFileNames = cms.untracked.vstring('file:/afs/cern.ch/user/k/koschwei/work/public/MuonEG_Run299368_v1_Run2017C_RAW_LS-79to90.root')
#)



process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
	#'/store/data/Run2017D/MuonEG/MINIAOD/17Nov2017-v1/70000/FAB825A5-ABE2-E711-845D-008CFA1C6448.root'
        'file:/afs/cern.ch/work/j/johnda/public/BJetTrigger/CMSSW_9_4_6_patch1/src/HLTBTagging/nTuples/crab/Run2017B_MuonEG_MINIAOD_9C2DFD61-5DE5-E711-BC68-20CF300E9ED0.root',
        #'/store/data/Run2017C/MuonEG/MINIAOD/17Nov2017-v1/60000/B20555F6-65E4-E711-B663-001CC4A6ABA8.root'
        #'file:/afs/cern.ch/user/k/koschwei/work/public/MuonEGRunC_MiniAOD_300107_3E580A66-3477-E711-8027-02163E0142F6.root'
    ),
    #lumisToProcess = cms.untracked.VLuminosityBlockRange("302031:1-302031:1"),
    #secondaryFileNames = cms.untracked.vstring('/store/data/Run2017D/MuonEG/RAW/v1/000/302/031/00000/74EA91B8-6C8D-E711-9382-02163E013521.root')
    lumisToProcess = cms.untracked.VLuminosityBlockRange("297047:4-297047:13"),
    secondaryFileNames = cms.untracked.vstring('file:/afs/cern.ch/work/j/johnda/public/BJetTrigger/CMSSW_9_4_6_patch1/src/HLTBTagging/nTuples/crab/Run2017B_MuonEG_RAW_86FB4930-DC52-E711-B81B-02163E01375F.root')
)

"""
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/mc/RunIIFall17MiniAOD/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/TSG_94X_mc2017_realistic_v11-v1/30000/AEA25818-FB21-E811-BBFB-848F69FD444B.root'),
    #lumisToProcess = cms.untracked.VLuminosityBlockRange("300107:2-300107:3"),
    lumisToProcess = cms.untracked.VLuminosityBlockRange("1:261-1:262"),
    secondaryFileNames = cms.untracked.vstring('/store/mc/RunIIFall17DRPremix/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/GEN-SIM-RAW/TSG_94X_mc2017_realistic_v11-v1/30000/14E63E79-041E-E811-A72A-3417EBE5289A.root')
)
"""
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

process.options = cms.PSet(
    numberOfThreads = cms.untracked.uint32(4)
)

process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('tree.root')
)


process.out = cms.EndPath(process.output)

