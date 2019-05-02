#!/usr/bin/env python
import os
#import PhysicsTools.HeppyCore.framework.config as cfg
#cfg.Analyzer.nosubdir=True

import sys
import re

import shlex
from subprocess import Popen, PIPE
def launch(cmd):
    process = Popen(shlex.split(cmd), stdout=PIPE)
    output = process.communicate() #output
    exit_code = process.wait()
    return exit_code

print "ARGV:",sys.argv
JobNumber=sys.argv[1]
islocal = False
try:
    ##CRAB
    import PSet
    crabFiles=PSet.process.source.fileNames
    crabSecondaryFiles=PSet.process.source.secondaryFileNames
    maxEvents=int(PSet.process.maxEvents.input.value())
    VLuminosityBlockRange = PSet.process.source.lumisToProcess
except:
    ##local test
    print "=================== I'm using local parameters ==================="
    import PSet_localTest as PSet
    crabFiles=PSet.process.source.fileNames
    crabSecondaryFiles=PSet.process.source.secondaryFileNames
    maxEvents=int(PSet.process.maxEvents.input.value())
    VLuminosityBlockRange = PSet.process.source.lumisToProcess
    islocal = True
if maxEvents<0:
    maxEvents = 1000000000
print "crabFiles before: ",crabFiles
print "crabSecondaryFiles before: ",crabSecondaryFiles

if len(crabFiles)>0:
    firstInput = crabFiles[0]
else:
    firstInput = "emptyFile"
print "crabFiles after: ",crabFiles
print "crabSecondaryFiles after: ",crabSecondaryFiles


from fwlite_config_phaseI import *
import imp
handle = open("fwlite_config_phaseI.py", 'r')
cfo = imp.load_source("config", "fwlite_config_phaseI.py", handle)
#config = cfo.config
handle.close()

##replace files with crab ones
#config.components[0].files=crabFiles

launchNtupleFromHLT("tree.root",crabFiles,crabSecondaryFiles,maxEvents, LS = VLuminosityBlockRange, local = islocal, preProcessing=True)

#from PhysicsTools.HeppyCore.framework.looper import Looper
#looper = Looper( 'Output', config, nPrint = 1)
#looper.loop()
#looper.write()

#print PSet.process.output.fileName
#os.system("ls -lR")
try:
    os.rename("Output/tree.root", "tree.root")
except:
    pass
#os.system("ls -lR")

import ROOT
try:
    f=ROOT.TFile.Open('tree.root')
    entries=f.Get('tree').GetEntries()
except:
    entries=0

fwkreport='''<FrameworkJobReport>
<ReadBranches>
</ReadBranches>
<PerformanceReport>
  <PerformanceSummary Metric="StorageStatistics">
    <Metric Name="Parameter-untracked-bool-enabled" Value="true"/>
    <Metric Name="Parameter-untracked-bool-stats" Value="true"/>
    <Metric Name="Parameter-untracked-string-cacheHint" Value="application-only"/>
    <Metric Name="Parameter-untracked-string-readHint" Value="auto-detect"/>
    <Metric Name="ROOT-tfile-read-totalMegabytes" Value="0"/>
    <Metric Name="ROOT-tfile-write-totalMegabytes" Value="0"/>
  </PerformanceSummary>
</PerformanceReport>
<GeneratorInfo>
</GeneratorInfo>
<InputFile>
<LFN>%s</LFN>
<PFN></PFN>
<Catalog></Catalog>
<InputType>primaryFiles</InputType>
<ModuleLabel>source</ModuleLabel>
<GUID></GUID>
<InputSourceClass>PoolSource</InputSourceClass>
<EventsRead>1</EventsRead>
</InputFile>
<File>
<LFN></LFN>
<PFN>tree.root</PFN>
<Catalog></Catalog>
<ModuleLabel>HEPPY</ModuleLabel>
<GUID></GUID>
<OutputModuleClass>PoolOutputModule</OutputModuleClass>
<TotalEvents>%s</TotalEvents>
<BranchHash>dc90308e392b2fa1e0eff46acbfa24bc</BranchHash>
</File>
</FrameworkJobReport>''' % (firstInput,entries)

f1=open('./FrameworkJobReport.xml', 'w+')
f1.write(fwkreport)
