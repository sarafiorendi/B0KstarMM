#####################################################################
# Program to fit with different parameter starting values:          #
# 1. generate new (random-)parameter files --> "name_bin_indx.txt"  #
# 2. use files at 1. to perform the fits (parFileName = "name.txt") #
#####################################################################

from os import system
import sys

if len(sys.argv) < 6:
    print "Synopsis: python runBatchCheckMuliMinima.py [FitType] [File.root] [FitEff] [q^2 bin] [nFiles] [For fits[runJobs[true / false]] [parFileName]]"
    sys.exit()

FitType  = sys.argv[1]
FileName = sys.argv[2]
FitEff   = sys.argv[3]
q2Bin    = sys.argv[4]
nFiles   = int(sys.argv[5])
runJobs  = "false"

if not (int(sys.argv[1]) == 96) and len(sys.argv) != 8:
    print "Synopsis: python runBatchCheckMuliMinima.py [FitType] [File.root] [FitEff] [q^2 bin] [nFiles] [For fits[runJobs[true / false]] [parFileName]]"
    sys.exit()
elif not(int(sys.argv[1]) == 96):
    runJobs = sys.argv[6]
    tmpName = sys.argv[7]

system("unset DISPLAY")

for i in xrange(0,nFiles):
    if int(sys.argv[1]) == 96:
        execProg = ".././ExtractYield " + FitType + " " + FileName + " " + FitEff + " " + q2Bin + " " + str(i)
    else:
        ParameterFILE = tmpName.replace(".txt","_" + q2Bin + "_" + str(i) + ".txt")
        execProg = ".././ExtractYield " + FitType + " " + FileName + " " + FitEff + " " + q2Bin + " " + str(i) + " " + ParameterFILE

    if runJobs == "true":
        system("sleep 2")
        toRun = "Qsub -l lnxfarm -e -o Fit_" + q2Bin + "_" + str(i) + ".log -N FIT" + q2Bin + str(i) + " " + execProg
    else:
        toRun = execProg

    print toRun
    system(toRun)
