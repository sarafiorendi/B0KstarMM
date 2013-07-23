###########################################################
# Program to fit with different parameter starting values #
###########################################################

from os import system
import sys

if len(sys.argv) < 6:
    print "Synopsis: python RunBatchCheckMuliMinima.py [FitType] [File.root] [FitEff] [q^2 bin] [nFiles] [[runJobs[true / false]] [parFileName]]"
    sys.exit()

FitType  = sys.argv[1]
FileName = sys.argv[2]
FitEff   = sys.argv[3]
q2Bin    = sys.argv[4]
nFiles   = int(sys.argv[5])
runJobs  = "false"

if not (int(sys.argv[1]) >= 93 and int(sys.argv[1]) <= 96) and len(sys.argv) != 8:
    print "Synopsis: python RunBatchCheckMuliMinima.py [FitType] [File.root] [FitEff] [q^2 bin] [nFiles] [[runJobs[true / false]] [parFileName]]"
    sys.exit()
elif not (int(sys.argv[1]) >= 93 and int(sys.argv[1]) <= 96):
    runJobs = sys.argv[6]
    tmpName = sys.argv[7]

system("unset DISPLAY")

for i in xrange(1,nFiles+1):
    if int(sys.argv[1]) >= 93 and int(sys.argv[1]) <= 96:
        execProg = "./ExtractYield " + FitType + " " + FileName + " " + FitEff + " " + q2Bin + " " + str(i)
    else:
        ParameterFILE = tmpName.replace(".txt","_" + q2Bin + "_" + str(i) + ".txt")
        execProg = "./ExtractYield " + FitType + " " + FileName + " " + FitEff + " " + q2Bin + " " + ParameterFILE + " " + str(i)

    if runJobs == "true":
        system("sleep 2")
        toRun = "Qsub -l lnxfarm -e -o Fit_" + q2Bin + "_" + str(i) + ".log -N FIT" + q2Bin + str(i) + " " + execProg
    else:
        toRun = execProg

    print toRun
    system(toRun)
